import argparse
import datetime
import multiprocessing
import os
import time  # for simple profiling
from decimal import *

# https://biopython.org/ for parsing structural information
from Bio.PDB import FastMMCIFParser, MMCIF2Dict

# https://www.sqlalchemy.org/ for generating SQL database
from sqlalchemy import create_engine
from sqlalchemy import exc
from sqlalchemy import literal
from sqlalchemy.orm import sessionmaker
from SQLalchemy_declarative import Base, DCDistance, DCSummary


def main():
    """
    This main function interprets the task and assigns analysis processes to each core.
    Data is written to a database via SQLalchemy.
    No return.
    """

    task = Task(args)  # generates a task object containing job information
    task.define_queue()  # defines the total queue of .cif PDB structure files to be analysed

    # Assign multiprocessing. Structure from https://sebastianraschka.com/Articles/2014_multiprocessing.html
    output = multiprocessing.Queue()
    processes = [multiprocessing.Process(target=pdb_process, args=(task, core, output)) for core in
                 range(task.cores)]  # split total file queue and assign processes with pdb_process function to cores
    for p in processes:
        p.start()
    results = [output.get() for p in processes]
    for result in results:
        print(result)
    for p in processes:
        p.join()
    print(processes)


class Task:  # Class managing information about the run itself, e.g. queue, parameters, arguments
    def __init__(self, args):  # loads arguments from textfile
        job_filepath = args.job  # filepath of textfile with arguments, defined by argparse args.job

        # Interprets input file. Structure from https://stackoverflow.com/questions/1305532/convert-nested-python-dict-to-object
        with open(job_filepath, "r") as job_fh:
            for row in job_fh:
                arg, value = row.replace("\n", "").split("\t")
                if arg == "cutoff":
                    value = (float(value))
                elif arg in ["chunksize", "chunksize_offset", "flush_offset", "modellimit", "cores", "filesize_limit"]:
                    value = int(value)
                elif arg in ["target_residue", "target_atom"]:
                    value = value.upper().replace(" ", "").replace('"', "")
                    value = value.split(",")
                if arg in ["target_NT", "distance_CT", "distance_NT", "permit_chain_identity", "verbose", "quiet",
                           "appending"]:
                    value = (value == "TRUE")
                self.__dict__[arg] = value

    def define_queue(self):
        if self.appending:
            return self.define_partial_queue()

        file_list = []
        for root, dirs, files in os.walk(self.pdb_database):
            for file in files:
                # Bugfix - Creating the pdb object here makes garbage collection impossible and RAM use explodes. Create during task only
                file_list.append((root, file))

        self.file_queue = []
        for i in range(0, self.cores):
            self.file_queue.append(file_list[i:len(file_list):self.cores])

        print("File queue of {} files split across {} cores.".format(len(file_list), self.cores))
        return True

    # This code allows for analysis of new structures without reanalysing the entire local PDB copy
    def define_partial_queue(self):  # TODO implement modified timestamp
        engine = create_engine(self.dc_db)
        Base.metadata.bind = engine
        db_session = sessionmaker(bind=engine)
        session = db_session()

        file_list = []
        dc_pass_count = 0
        with open(self.dc_pass_log, "w+") as outfile:
            for root, dirs, files in os.walk(self.pdb_database):
                for file in files:
                    #  Bugfix - Creating the pdb object here makes garbage collection impossible and RAM use explodes. Create during task only
                    file_id = file.replace(".cif", "")
                    q = session.query(DCSummary).filter(DCSummary.pdb_id == file_id)
                    # print(session.query(q.exists()))
                    if (session.query(literal(True)).filter(q.exists()).scalar()):
                        outfile.write(file_id + "\n")
                        dc_pass_count += 1
                    else:
                        file_list.append((root, file))

        self.file_queue = []
        for i in range(0, self.cores):
            self.file_queue.append(file_list[i:len(file_list):self.cores])

        print("Partial file queue of {} files, skipped {} files. Split across {} cores.".format(len(file_list),
                                                                                                dc_pass_count,
                                                                                                self.cores))
        session.close()
        return True


# This is the target function for multiprocessing which handles analysis of the structure file queue
def pdb_process(task, core, output):
    process_start_time = time.time()  # for simple profiling

    chunksize = task.chunksize + task.chunksize_offset * core  # chunksize: how many files to analyse before commiting data to the database

    # assigning SQLalchemy connection
    engine = create_engine(task.dc_db)
    Base.metadata.bind = engine
    db_session = sessionmaker(bind=engine)
    session = db_session()
    str_parser = FastMMCIFParser(QUIET=1)

    print("Connected core %i, chunksize %i, offset %i" % (core, task.chunksize, task.chunksize_offset * core))

    # detail of console output
    verbose = task.verbose
    quiet = task.quiet

    # buffer to hold data for submission to database and if it has been submitted successfully. Handles potential of locked database while writing from parallel process
    buffer = []
    flush = True
    flush_offset = task.flush_offset
    flush_offset_count = 0
    file_queue = task.file_queue[core]

    # loop to analyse individual structure files from file queue
    for pdb_count, (root, file) in enumerate(file_queue):
        pdb = PDB(os.path.join(root, file),
                  core)  # generate a PDB object to hold all relevant information for a single parsed structure file
        task_summary = ""
        task_summary = task_summary + ("Processing file %s on core %i" % (pdb.path, core))

        # checks before structure analysis, including successful fast structure parsing
        if pdb.filesize > task.filesize_limit:  # throw out large files. TODO: Seperate filesize due to large assembly from structure factors
            task_summary = task_summary + "\n\tSize abort"
            pdb.pass_filesize = False
            pdb.abort = True
            print(task_summary)
        else:
            if task.verbose: task_summary = task_summary + "\n\tSize pass"
            pdb.pass_filesize = True
            try:
                pdb.structure = str_parser.get_structure(pdb.id, pdb.path)
                if task.verbose: task_summary = task_summary + "\n\tMMCIFParser pass"
                pdb.pass_structure = True
            except:
                task_summary = task_summary + "\n\tMMCIFParser abort"
                pdb.pass_structure = False
                pdb.abort = True
                print(task_summary)

        # PDB object functions analysing distance data in a structure file and assigning relevant metadata for the database
        pdb.analyse(task)
        pdb.generate_dictionary()
        pdb.assign_dictionary_data()
        pdb.sqla_convert_distances()

        if len(
                pdb.filtered_distances) > 0:  # hits from filtered_distances are preferred the summary table (here: intermolecular hits)
            representative_distance = pdb.representative_distance_filtered()
        else:  # no intermolecular hits found or no hits found at all
            representative_distance = pdb.representative_distance_unfiltered()
        if representative_distance:
            pdb.top_hit_sqla(representative_distance)  # populate class with information for summary table

        pdb.sqla_summary()  # generate entry for summary table

        pdb_entries_wrap = pdb.alchemy_distances, pdb.alchemy_sum  # results for all distances per structure and single distance for summary table

        buffer.append(pdb_entries_wrap)  # append result to buffer holding data for submission to database
        if len(buffer) == chunksize:
            print("Core %i - %i tasks in %i seconds" % (core, pdb_count + 1, (time.time() - process_start_time)))
            flush = buffer_SQLal_dc_submission(session, buffer, core)  # returns if submission to database successful
            chunksize = task.chunksize  # remove offset after first time. This might not be needed

            if flush:
                flush_offset_count = 0
                buffer = []
                print("Core %i - commit successful" % core)

        if not flush:  # Handles incomplete submission of buffer due to busy database. This fix likely already made chunk_offset obsolete
            if flush_offset_count == flush_offset:  # Tries to submit again every flush_offset until it succeeds
                print("Core %i - not flush at %i tasks with %i tasks in buffer" % (core, pdb_count + 1, len(buffer)))
                flush_offset_count = 0
                flush = buffer_SQLal_dc_submission(session, buffer, core)
                if flush:
                    buffer = []
            flush_offset_count += 1

        # # # some text for console
        if not verbose:
            task_summary = task_summary + "\nFinished task %i on core %i\n" % (pdb_count, core)
        if verbose:
            task_summary = (task_summary + ("\n\t\t%i\tdistances\t" % len(pdb.distances) + str(pdb.distances)))
            task_summary = task_summary + "\n\t\t\t%i of %i distances intermolecular" % (
                pdb.inter_count, len(pdb.distances))
            task_summary = task_summary + "\n\t\t\t%i of %i distances below cutoff" % (
                pdb.dist_count, len(pdb.distances))
            task_summary = task_summary + "\n\t\t\t%i distances below cutoff and intermolecular" % (pdb.hit_count)
            task_summary = task_summary + "\nFinished task %i on core %i\n" % (pdb_count, core)
        if not quiet: print(task_summary)

    #  Final commit for incomplete buffer at end of queue (remainder from chunk size)
    for i in range(0, 10):
        flush = buffer_SQLal_dc_submission(session, buffer, core)  # final commit attempted 10 times
        if not flush:
            time.sleep(60)
            print("Failed to commit final chunk of size %i on core %i. Waiting 60s" % (len(buffer), core))
        else:
            print("Committed final chunk of size %i on core %i." % (len(buffer), core))
            break
    if not flush:
        print("Failed to flush after 10 min. Attempting single submission")
        for i in range(0, 10):
            buffer = single_buffer_SQLal_dc_submission(session, buffer, core)  # final commit attempted 10 times
            if len(buffer) == 0:
                flush = True
            else:
                flush = False

            if not flush:
                time.sleep(60)
                print("Failed to single commit, final chunk of size %i on core %i. Waiting 60s" % (len(buffer), core))
            else:
                print("Committed all chunks on core %i." % (core))
                break
    if not flush:  # writes error to logfile upon failure to commit final entries after 2*10 attempts
        with open(task.dc_error_log, "a+") as outfile:
            for entry in buffer:
                outfile.write("{}\t{}".format(entry[1].pdb_id, core))

    output.put("--- Core %s finished %s tasks in %s seconds ---" % (
        core, len(file_queue), (time.time() - process_start_time)))


def buffer_SQLal_dc_submission(session, buffer, core):  # commits entries in buffer to database
    """
    Args:
        session: SQLalchemy session connecting to database
        buffer: List containing pdb_wrap tuples
    Returns:
        bool: submission to database successful
    """

    # # Add entries to session
    for pdb_wrap in buffer:
        for alchemy_entry in pdb_wrap[0]:  # [0]: pdb.alchemy_entries
            session.add(alchemy_entry)
        session.add(pdb_wrap[1])

    # # Try to submit entries to database
    try:
        #  session.flush()
        session.commit()
        return True

        # # Failure to submit when database is locked, e.g. busy (other threads submitting)
    except exc.OperationalError as e:
        #  Rolls back database to try again later (see pdb_process)
        print("SQLalchemy exc.OperationalError on core %i, buffer size %i. Rolling back database. Error message: %s" % (
            core, len(buffer), e))
        session.rollback()
        return False


def single_buffer_SQLal_dc_submission(session, buffer, core):  # commits entries in buffer to database, entry by entry
    # # Add entries to session
    buffer2 = []
    for pdb_wrap in buffer:
        for alchemy_entry in pdb_wrap[0]:  # [0]: pdb.alchemy_entries
            session.add(alchemy_entry)
        session.add(pdb_wrap[1])
        try:
            #  session.flush()
            session.commit()
        except exc.OperationalError as e:
            #  Rolls back database to try again later (see pdb_process)
            print(
                "Single submission SQLalchemy exc.OperationalError on core %i, buffer size %i. Rolling back database. Full error %s" % (
                    core, len(pdb_wrap), e))
            session.rollback()
            buffer2.append(pdb_wrap)
    return buffer2


class PDB:  # main class to handle structure interpretation, analysis, metadata
    #  def __del__(self):
    #     print("PDB deleted")

    def __init__(self, filepath, core):
        self.core = core

        #  Functionality
        self.path = filepath
        self.filesize = os.stat(filepath).st_size
        self.id = os.path.split(filepath)[1].replace(".cif", "")
        if "-" in filepath:  # determines if pdb file was constructed from bio assembly code or not
            self.bio = int(filepath.split("-")[1].split(".")[0])
        else:
            self.bio = None

        self.structure = []
        self.filtered_distances = []
        self.distances = []

        #  Set SQL summary entry defaults
        self.pass_filesize = None
        self.pass_structure = None
        self.pass_dictionary = None
        self.pass_analysis = None
        self.hit_analysis = None
        self.abort = False
        self.dist = None
        self.cys_count = 0  # count of cysteines
        self.dist_count = 0  # count of distances below cutoff
        self.inter_count = 0  # count of intermolecular distances
        self.hit_count = 0  # count of intermolecular distances below cutoff

        #  Set SQLalchemy distance object defaults
        self.nt = None
        self.d1_X_C = None
        self.d2_X_CA = None
        self.d3_X_N = None
        self.d4_CA_C = None
        self.d5_CA_CA = None
        self.d6_CA_N = None
        self.d7_N_C = None
        self.d8_N_CA = None
        self.d9_N_N = None
        self.x_name = None
        self.dynamic_chain = None
        self.dynamic_resn = None
        self.dynamic_resi = None
        self.dynamic_pos = None
        self.dynamic_name = None
        self.dynamic_order = None
        self.static_chain = None
        self.static_resn = None
        self.static_resi = None
        self.static_pos = None
        self.static_name = None
        self.static_order = None
        self.model_count = None
        self.static_chain_res_count = None
        self.poly_ids = None
        self.dynamic_chain_res_count = None
        self.model = None
        self.chain_count = None
        self.distances_count = None  # count of all calculated distances
        self.dynamic_poly_id = None
        self.static_poly_id = None
        self.poly_id_eq = None
        self.inter = None
        self.timestamp = str(datetime.datetime.utcnow())

    def analyse(self, task):  # calls functions to analyse distances for given structure
        self.job_id = task.jobname
        self.model_count = len(self.structure)

        if self.model_count > 0:
            count = 0
            for model in self.structure:
                if count < task.modellimit:
                    self.define_residues(task, model)
                    self.generate_distances(task)
                else:
                    break
                count += 1
            self.evaluate_distances(task)
            self.average_models(count)
            self.pass_analysis = True

    def define_residues(self, task,
                        model):  # identifies 'static' (e.g. C-terminus) and 'dynamic'/target (e.g. Lysines) residues
        target_list, nts, cts = [], [], []
        chain_count = 0
        for chain in model:
            chain.index = chain_count
            res_count = 0  # counts resolved residues in chain
            ct, nt = None, False
            for pos, residue in enumerate(chain):
                if not residue.get_resname() in aminoacid: continue  # skip heteroresidues
                if not (residue.id[0] == ' '): continue  # skip heteroresidues

                res_count += 1
                ct = ((pos, residue))  # assigns ct to current residue until end of list, skipping heteroresidues
                if nt is False:
                    nt = True
                    nts.append((pos, residue, True))  # True/False: is N-terminus?

                #  can cause duplicate errors if kept together (N-terminus = Lysine if no atoms resolved in either)
                else:
                    if residue.get_resname() in task.target_residue:
                        target_list.append((pos, residue, False))
                if residue.get_resname() == 'CYS':  # tracks cysteines as these can be of special interest for protein production
                    self.cys_count += 1

            chain.res_count = res_count
            if res_count > 0:
                chain_count += 1
            if ct is not None:
                cts.append(ct)
        model.chain_count = chain_count
        self.dynamic_list = target_list + nts  # resets everytime
        self.static_list = cts  # resets everytime

    def average_models(self, model_count):  # determines average values for analysis of multiple models
        self.distances_count = Decimal("%.3f" % (self.distances_count / model_count))
        self.dist_count = Decimal("%.3f" % (self.dist_count / model_count))
        self.inter_count = Decimal("%.3f" % (self.inter_count / model_count))
        self.hit_count = Decimal("%.3f" % (self.hit_count / model_count))
        self.cys_count = Decimal("%.3f" % (self.cys_count / model_count))

    def generate_distances(self, task):  # calculates distances from residues in static list to residues in dynamic list
        for static_pos, static_res in self.static_list:
            for dynamic_pos, dynamic_res, nt_flag in self.dynamic_list:
                distance = Distance(self, dynamic_pos, dynamic_res, task.target_atom, nt_flag, static_pos, static_res)
                self.distances.append(distance)
        self.distances_count = len(self.distances)

    def evaluate_distances(self,
                           task):  # Pre-evaluates distances for summary table. These are partially unnecessary given subsequent SQL analysis of the distance table
        for distance in self.distances:
            distance.assign_cutoff(task)
            distance.assign_intra()
            distance.assign_hit()
            if not distance.intra:  # More conservatively this could also be distance.hit (considering cutoff)
                self.filtered_distances.append(distance)

        for distance in self.distances:
            if distance.cut:
                self.dist_count += 1  # only distances below cutoff, distances_count tracks all distances
                if not distance.intra:
                    self.hit_count += 1
            if not distance.intra:
                self.inter_count += 1

    def sqla_convert_distances(self):
        self.alchemy_distances = []
        for distance in self.distances:
            distance.assign_chain_equality()  # requires dictionary to be converted first
            distance.sqlalchemy_conversion()
            self.alchemy_distances.append(distance.alchemy)

    def sort_distances(self):
        self.distances = sorted(self.distances)

    def sort_filtered_distances(self):
        self.filtered_distances = sorted(self.filtered_distances)

    def generate_dictionary(self):
        self.pdb_dict = MMCIF2Dict.MMCIF2Dict(self.path)
        self.name = self.convert_dictionary_data('_struct.title')
        if self.name is None:
            self.pass_dictionary = False
        else:
            self.pass_dictionary = True

    def assign_dictionary_data(self):
        self.bio_list = self.convert_dictionary_data('_pdbx_struct_assembly_gen.asym_id_list')
        self.head = self.convert_dictionary_data('_struct.pdbx_descriptor')
        self.genus = self.convert_dictionary_data('_entity_src_gen.pdbx_gene_src_scientific_name')
        self.genes = self.convert_dictionary_data('_entity_src_gen.pdbx_gene_src_gene')
        self.host = self.convert_dictionary_data('_entity_src_gen.pdbx_host_org_scientific_name')
        self.deposition = self.convert_dictionary_data('_pdbx_database_status.recvd_initial_deposition_date')
        self.bio_count = self.count_bio()
        self.method = self.convert_dictionary_data('_exptl.method')
        self.resolution = self.convert_dictionary_data('_refine_hist.d_res_high')

        self.poly_ids = self.convert_dictionary_data('_entity_poly.pdbx_strand_id')
        if self.poly_ids is not None:
            #  converts string to array
            self.poly_ids = self.poly_ids.replace("[", "").replace("]", "").replace(" ", "").split("'")
            idsl = []
            for ids in self.poly_ids:
                ids = ids.split(",")
                ids_e = []
                for element in ids:
                    if len(element) > 0:
                        ids_e.append(element)
                if len(ids_e) > 0:
                    idsl.append(ids_e)
            self.poly_ids = idsl

    def convert_dictionary_data(self, key):
        entry = str(self.pdb_dict.get(key))
        if entry == "None":
            entry = None
        return entry

    def count_bio(self):  # Determines average count of chains in biological assembly
        #  This is only somewhat informative as it includes 'chains' which are just heteroatoms / water /...
        entry = self.pdb_dict.get('_pdbx_struct_assembly_gen.asym_id_list')
        if type(entry) == list:
            i = 0
            for j in entry:
                i = i + len(j.replace(',', ""))
            i = i / len(entry)
            return Decimal("%.3f" % i)
        elif isinstance(entry, str):
            return Decimal("%.3f" % len(entry.replace(',', "")))
        return None

    def representative_distance_filtered(
            self):  # determines a representative distance to use to populate summary table. Herein prioritizing intermolecular
        self.sort_filtered_distances()
        #  search for first distance between non-bio chains and of priority nucleophile
        index = None
        alt1index = None
        alt2index = None

        #  a bit complex way to find "min" -> Might be easier to split into seperate tables by condition instead, then determine min for each
        for pos, distance in enumerate(
                self.filtered_distances):  # identifies first filtered distance with non-equal chains
            if not distance.poly_id_eq and distance.poly_id_eq is not None:  # or (distance.ref_id_eq is False): /ref_id is not reliable for biological assemblies due to reassignment of chains
                if (distance.dynamic_res.get_resname() == "LYS") or distance.nt:
                    index = pos
                    break  # break at first priority nucleophile satisfying not distance.poly_id_eq
                elif alt1index is None:  # first one satisfying not distance.poly_id_eq which is not also priority nucleophile
                    altindex = pos
                    continue
            elif alt2index is None:  # first one satisfying condition
                if (distance.dynamic_res.get_resname() == "LYS") or distance.nt:
                    alt2index = pos

        if index is None:
            if alt1index is not None:
                index = alt1index
            else:
                if alt2index is not None:
                    index = alt2index
                else:
                    index = 0

        return self.filtered_distances[index]  # top_distance

    def representative_distance_unfiltered(self):
        dis_list = []
        for distance in self.distances:
            dis = distance.distance
            if dis is not None:
                dis_list.append(dis)
        if len(dis_list) == 0:
            return False
            #  https://stackoverflow.com/questions/2474015/getting-the-index-of-the-returned-max-or-min-item-using-max-min-on-a-_list
        index = min(range(len(dis_list)), key=dis_list.__getitem__)  # Finds index for minimum distance
        return self.distances[index]  # top_distance

    def top_hit_sqla(self, top_distance):  # populates information based on top distance
        self.dist = top_distance.distance
        self.model = top_distance.static_res.parent.parent.id
        self.chain_count = top_distance.static_res.parent.parent.chain_count
        self.dynamic_resn = str(top_distance.dynamic_res.get_resname())
        self.static_resn = str(top_distance.static_res.get_resname())
        self.dynamic_chain = str(top_distance.dynamic_res.parent.id)
        self.static_chain = str(top_distance.static_res.parent.id)
        self.dynamic_resi = top_distance.dynamic_res.id[1]
        self.dynamic_pos = top_distance.dynamic_pos
        self.static_resi = top_distance.static_res.id[1]
        self.dynamic_order = top_distance.dynamic_res.disordered
        self.static_order = top_distance.static_res.disordered
        self.static_pos = top_distance.static_pos
        self.static_chain_res_count = top_distance.static_res.parent.res_count
        self.dynamic_chain_res_count = top_distance.dynamic_res.parent.res_count
        self.dynamic_poly_id = str(top_distance.dynamic_res.parent.poly_id)
        self.static_poly_id = str(top_distance.static_res.parent.poly_id)
        self.poly_id_eq = top_distance.poly_id_eq
        self.nt = top_distance.nt
        self.d1_X_C = top_distance.d1_X_C
        self.d2_X_CA = top_distance.d2_X_CA
        self.d3_X_N = top_distance.d3_X_N
        self.d4_CA_C = top_distance.d4_CA_C
        self.d5_CA_CA = top_distance.d5_CA_CA
        self.d6_CA_N = top_distance.d6_CA_N
        self.d7_N_C = top_distance.d7_N_C
        self.d8_N_CA = top_distance.d8_N_CA
        self.d9_N_N = top_distance.d9_N_N
        self.x_name = top_distance.x_name
        self.inter = not top_distance.intra

        if top_distance.dynamic is None:
            print("filtered_distances[index].dynamic is None on file %s on core %i" % (self.id, self.core))
            #  Note this can occur for unresolved N-terminus
            self.dynamic_name = None
        else:
            self.dynamic_name = str(top_distance.dynamic.get_name())

        if top_distance.static is None:
            self.static_name = None
            #  Note this can occur for bad file, e.g. 3nso CYS'99 only SG resolved
        else:
            self.static_name = str(top_distance.static.get_name())
        #  if args.verbose: print("hit_count", len(self.alchemy_entries), "dist_count", self.dist_count, "filtered_distances", len(self.filtered_distances))

    def sqla_summary(self):
        self.alchemy_sum = DCSummary(
            abort=self.abort,
            pass_filesize=self.pass_filesize,
            filesize=self.filesize,
            pass_structure=self.pass_structure,
            pass_dictionary=self.pass_dictionary,
            pass_analysis=self.pass_analysis,
            dynamic_chain=self.dynamic_chain,
            dynamic_chain_res_count=self.dynamic_chain_res_count,
            dynamic_resn=self.dynamic_resn,
            dynamic_resi=self.dynamic_resi,
            dynamic_name=self.dynamic_name,
            dynamic_order=self.dynamic_order,
            static_chain=self.static_chain,
            static_chain_res_count=self.static_chain_res_count,
            static_resn=self.static_resn,
            static_resi=self.static_resi,
            static_pos=self.static_pos,
            static_name=self.static_name,
            static_order=self.static_order,
            model_count=self.model_count,
            model=self.model,
            name=self.name,
            head=self.head,
            genes=self.genes,
            genus=self.genus,
            host=self.host,
            deposition=self.deposition,
            method=self.method,
            resolution=self.resolution,
            bio=self.bio,
            bio_list=self.bio_list,
            bio_count=self.bio_count,
            chain_count=self.chain_count,
            pdb_id=self.id,
            distance=self.dist,
            d1_X_C=self.d1_X_C,
            d2_X_CA=self.d2_X_CA,
            d3_X_N=self.d3_X_N,
            d4_CA_C=self.d4_CA_C,
            d5_CA_CA=self.d5_CA_CA,
            d6_CA_N=self.d6_CA_N,
            d7_N_C=self.d7_N_C,
            d8_N_CA=self.d8_N_CA,
            d9_N_N=self.d9_N_N,
            nterm=self.nt,
            x_name=self.x_name,
            distances=self.distances_count,
            poly_ids=str(self.poly_ids),
            dynamic_chain_poly_id=self.dynamic_poly_id,
            static_chain_poly_id=self.static_poly_id,
            chain_poly_id_eq=self.poly_id_eq,
            inter=self.inter,
            dist_count=self.dist_count,
            inter_count=self.inter_count,
            hit_count=self.hit_count,
            cys_count=self.cys_count,
            job_id=self.job_id,
            timestamp=self.timestamp
        )


class Distance:  # Class to handle distances between two residues
    #  consider give it function to test different atoms itself instead of PDB - res-distance instead of distance?
    #  def __del__(self):
    #     print("Distance deleted")

    def __init__(self, pdb, dynamic_pos, dynamic_res, dynamic_target_atoms, nt_flag, static_pos, static_res):
        self.parent_pdb = pdb
        self.static_pos = static_pos
        self.dynamic_pos = dynamic_pos
        self.dynamic_res = dynamic_res
        self.static_res = static_res

        #  Set SQLalchemy distance object defaults
        self.nt = nt_flag
        self.d1_X_C = None
        self.d2_X_CA = None
        self.d3_X_N = None
        self.d4_CA_C = None
        self.d5_CA_CA = None
        self.d6_CA_N = None
        self.d7_N_C = None
        self.d8_N_CA = None
        self.d9_N_N = None
        self.dynamic = None
        self.x_name = None

        #  fill values for d1 to d9
        self.calculate_distances(dynamic_target_atoms)
        #  Assign static atom to C, then CA, then N, then None
        self.assign_static()
        #  This code assigns the primary distance for conflicts between X and N in case of N-terminal residue (favouring N for no static atom)
        self.assign_dynamic()
        #  This determines the primary distance
        self.distance = self.atom_atom_distance(self.dynamic, self.static)

        #  adjust from array index to actual count
        self.static_pos += 1
        self.dynamic_pos += 1

    def calculate_distances(self, dynamic_target_atoms):
        for dynamic_atom in self.dynamic_res:
            if dynamic_atom.id in dynamic_target_atoms:  # Note: single nucleophile by residue / should be unique -> Lookup table by residue instead?
                self.x_name = dynamic_atom.id  # Note: this should be potentially be assigned before, not by searching in dynamic_atoms (if dynamic_atom.id == self.x_name)
                self.dynamic = dynamic_atom  # Assigns X as default for dynamic, can change in self.assign_dynamic()
                self.d1_X_C = self.atom_residue_distance(dynamic_atom, self.static_res, 'C')
                self.d2_X_CA = self.atom_residue_distance(dynamic_atom, self.static_res, 'CA')
                self.d3_X_N = self.atom_residue_distance(dynamic_atom, self.static_res, 'N')

            elif dynamic_atom.id == 'CA':
                self.d4_CA_C = self.atom_residue_distance(dynamic_atom, self.static_res, 'C')
                self.d5_CA_CA = self.atom_residue_distance(dynamic_atom, self.static_res, 'CA')
                self.d6_CA_N = self.atom_residue_distance(dynamic_atom, self.static_res, 'N')

            elif dynamic_atom.id == 'N':
                self.d7_N_C = self.atom_residue_distance(dynamic_atom, self.static_res, 'C')
                self.d8_N_CA = self.atom_residue_distance(dynamic_atom, self.static_res, 'CA')
                self.d9_N_N = self.atom_residue_distance(dynamic_atom, self.static_res, 'N')

    def assign_static(self):  # assigns preferred static atom (C > CA > N)

        for key in ('C', 'CA', 'N'):
            if key in self.static_res:
                self.static = self.static_res[key]
                break
        else:
            self.static = None
            self.static_name = None
            return False

        self.static_name = self.static.get_name()
        return True

    def assign_dynamic(self):  # assigns preferred dynamic atom (also see methods in corresponding paper)
        if self.nt:  # for N-terminal residue, re-evaluate 'dynamic' even if it already exists (can be pre-defined in class function calculate_distances)
            if self.static is None:
                try:
                    self.dynamic = self.dynamic_res[
                        'N']  # With no distances, default self.dynamic to N before X (best nucleophile)
                except KeyError:
                    pass
            else:
                if self.static.id == 'C':
                    if self.d7_N_C is None:  # If N doesn't exist, default to X
                        pass
                    elif self.d1_X_C is not None:  # If both N and X exist, compare distances and choose shorter one
                        if self.d1_X_C > self.d7_N_C:
                            self.dynamic = self.dynamic_res['N']
                    else:  # If X doesn't exist, default to N
                        self.dynamic = self.dynamic_res['N']
                elif self.static.id == 'CA':
                    if self.d8_N_CA is None:
                        pass
                    elif self.d2_X_CA is not None:
                        if self.d2_X_CA > self.d8_N_CA:
                            self.dynamic = self.dynamic_res['N']
                    else:  # If X doesn't exist, default to N
                        self.dynamic = self.dynamic_res['N']
                elif self.static.id == 'N':
                    if self.d9_N_N is None:
                        pass
                    elif self.d3_X_N is not None:
                        if self.d3_X_N > self.d9_N_N:
                            self.dynamic = self.dynamic_res['N']
                    else:  # If X doesn't exist, default to N
                        self.dynamic = self.dynamic_res['N']

        if self.dynamic is None:  # If X doesn't exist and dynamic wasn't found for N-terminus either, default to other positions

            if 'CA' in self.dynamic_res:
                self.dynamic = self.dynamic_res['CA']
                self.dynamic_name = self.dynamic.get_name()
                return True
            elif 'N' in self.dynamic_res:
                self.dynamic = self.dynamic_res['N']
                self.dynamic_name = self.dynamic.get_name()
                return True
            else:
                self.dynamic = None
                self.dynamic_name = None
                return False
        else:
            self.dynamic_name = self.dynamic.get_name()
            return True

    def atom_atom_distance(self, atom_1, atom_2):
        if atom_1 is not None and atom_2 is not None:
            return Decimal("%.3f" % (atom_1 - atom_2))
        else:
            return None

    def atom_residue_distance(self, atom_1, residue_2, target_2):
        """Determines a distance with None if atom not in residue"""
        if residue_2.__contains__(target_2):
            return Decimal("%.3f" % (atom_1 - residue_2[target_2]))
        else:
            return None

    def residue_residue_distance(self, residue_1, target_1, residue_2, target_2):
        """Determines a distance with None if atom not in residue"""
        if residue_1.__contains__(target_1):
            if residue_2.__contains__(target_2):
                return Decimal("%.3f" % (residue_1[target_1] - residue_2[target_2]))
            else:
                return None
        else:
            return None

    def assign_cutoff(self, task):
        if self.distance:
            if self.distance <= task.cutoff:
                self.cut = True
                return True
            else:
                self.cut = False
                return False
        else:
            self.cut = None  # return None

    def assign_intra(self):
        #  Chain identity
        if self.static_res.get_parent() is self.dynamic_res.get_parent():
            self.intra = True
            return True
        else:
            self.intra = False
            return False

    def assign_hit(self):
        if self.cut and not self.intra:
            self.hit = True
            return True
        else:
            self.hit = False
            return False

    def assign_chain_equality(self):  # requires dictionary information
        self.assign_chain_poly_id(self.get_parent_pdb().poly_ids)

        if self.dynamic_res.parent.poly_id == self.static_res.parent.poly_id:
            self.poly_id_eq = True
        else:
            self.poly_id_eq = False

    def assign_chain_poly_id(self, poly_ids):
        for residue in self.dynamic_res, self.static_res:
            try:
                residue.parent.poly_id
                #  print("residue.parent.poly_id already assigned")
            except AttributeError:
                if poly_ids is None:
                    residue.parent.poly_id = None
                else:
                    for ids_p, ids in enumerate(poly_ids):
                        if residue.parent.id in ids:
                            residue.parent.poly_id = ids_p

    def get_parent_pdb(self):
        return self.parent_pdb

    def sqlalchemy_conversion(self):
        self.alchemy = DCDistance(
            pdb_id=self.parent_pdb.id,
            parent_id=self.parent_pdb.id,
            distance=self.distance,
            d1_X_C=self.d1_X_C,
            d2_X_CA=self.d2_X_CA,
            d3_X_N=self.d3_X_N,
            d4_CA_C=self.d4_CA_C,
            d5_CA_CA=self.d5_CA_CA,
            d6_CA_N=self.d6_CA_N,
            d7_N_C=self.d7_N_C,
            d8_N_CA=self.d8_N_CA,
            d9_N_N=self.d9_N_N,
            nterm=self.nt,
            model=str(self.static_res.parent.parent.id),
            dynamic_resn=str(self.dynamic_res.get_resname()),
            dynamic_name=self.dynamic_name,
            x_name=self.x_name,
            static_name=self.static_name,
            static_resn=str(self.static_res.get_resname()),
            dynamic_chain=str(self.dynamic_res.parent.id),
            dynamic_chain_poly_id=str(self.dynamic_res.parent.poly_id),
            static_chain_poly_id=str(self.static_res.parent.poly_id),
            poly_id_eq=self.poly_id_eq,
            dynamic_chain_res_count=self.dynamic_res.parent.res_count,
            static_chain=str(self.static_res.parent.id),
            static_chain_res_count=self.static_res.parent.res_count,
            dynamic_resi=self.dynamic_res.id[1],
            dynamic_pos=self.dynamic_pos,
            static_resi=self.static_res.id[1],
            static_pos=self.static_pos,
            dynamic_order=self.dynamic_res.disordered,
            static_order=self.static_res.disordered,
            name=self.get_parent_pdb().name,
            inter=not self.intra,
            job_id=self.get_parent_pdb().job_id,
            timestamp=self.parent_pdb.timestamp  # almost no difference from self.timestamp = datetime.datetime.utcnow()
        )

    def __lt__(self, other):  # enables sorting of pdb.distances
        if self.distance is None:
            return False  # sorts None to highest position
        if other.distance is None:
            return True  # sorts None to highest position
        else:
            return self.distance < other.distance

    def __repr__(self):
        return "<dis %s, intra %s, hit %s>" % (self.distance, self.intra, self.hit)


aminoacid = ['VAL', 'ILE', 'LEU', 'GLU', 'GLN', 'ASP', 'ASN', 'HIS', 'TRP', 'PHE', 'TYR', 'ARG', 'LYS', 'SER', 'THR',
             'MET', 'ALA', 'GLY', 'PRO', 'CYS']  # defines allowed aminoacid names

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--job")
    parser.add_argument("-v", "--verbose", default=False)
    args = parser.parse_args()

    start_time = time.time()
    main()
    print("--- Program finished in %s seconds ---" % (time.time() - start_time))
