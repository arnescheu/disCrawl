import csv
import time

from sqlalchemy import create_engine
from sqlalchemy import exc
from sqlalchemy.orm import sessionmaker

from SQLalchemy_declarative import Base, DCSummary, CompDistanceUnique, Sub2DistanceUnique

def main(query_file, db_input, db_output,
         log_file_csv):  # query_file: tab seperated 'key|entry' e.g. 'pdb_id	|5glh-1.bio_str'
    queue = define_queue(query_file)
    subset_database(queue, db_input, db_output, log_file_csv)


def define_queue(input_file):
    queue = []
    with open(input_file, 'r') as csvfile:
        logreader = csv.reader(csvfile, delimiter="|")
        for row in logreader:
            queue.append((row[0].strip(), row[1].strip()))
    return queue


def subset_database(queue, db_input, db_output, log_file_csv):
    reset_database(db_output)

    db_in = create_engine(db_input)
    db_out = create_engine(db_output)
    DB_in_Session = sessionmaker(db_in)
    DB_out_Ssession = sessionmaker(db_out)
    in_session = DB_in_Session()
    out_session = DB_out_Ssession()

    with open(log_file_csv, "w+") as out_file:
        out_file.write("Count|Total|Key|Keyword|Instances|Submitted\n")
        for count, entry in enumerate(queue):
            print("Query {} of {} \t- {}".format(count, len(queue), entry))
            q = switch_key(entry, in_session)
            q_count = 0
            s_count = 0
            for instance in q:
                summary = rebuild_summary(instance, entry)
                if submit_summary(summary, out_session):
                    s_count += 1
                q_count += 1
            out_file.write("{}|{}|{}|{}|{}|{}\n".format(count, len(queue), entry[0], entry[1], q_count, s_count))
            out_file.flush()
            print("\tReceived {} instances, submitted {} new entries.".format(q_count, s_count))

    in_session.close()
    out_session.close()


def reset_database(db_output):
    print("Deleting {}".format(db_output.replace("sqlite:///", "")))
    # os.remove(db_output.replace("sqlite:///",""))
    engine = create_engine('sqlite:///C://Users//Arne//Desktop//disCrawl_filtered.db')
    Base.metadata.create_all(engine)
    return True


def switch_key(entry, in_session):
    key = entry[0]
    term = entry[1]
    if key == "pdb_id":
        q = in_session.query(DCSummary).filter(DCSummary.pdb_id.like("%{}%".format(term)))
    elif key == "name":
        q = in_session.query(DCSummary).filter(
            DCSummary.name.like("%{}%".format(term)))  # Add whitespaces around to avoid midword hits (e.g. "MET")
    elif key == "head":
        q = in_session.query(DCSummary).filter(
            DCSummary.head.like("%{}%".format(term)))  # Add whitespaces around to avoid midword hits (e.g. "MET")
    # TODO needs to throw error if undefined, or handle otherwise
    return q

def rebuild_summary(query, entry):  # TODO create seperate database and/or add query to filtered entry
    alchemy_sum = DCSummary(
        abort=query.abort,
        pass_filesize=query.pass_filesize,
        filesize=query.filesize,
        pass_structure=query.pass_structure,
        pass_dictionary=query.pass_dictionary,
        pass_analysis=query.pass_analysis,
        dynamic_chain=query.dynamic_chain,
        dynamic_chain_res_count=query.dynamic_chain_res_count,
        dynamic_resn=query.dynamic_resn,
        dynamic_resi=query.dynamic_resi,
        dynamic_name=query.dynamic_name,
        dynamic_order=query.dynamic_order,
        static_chain=query.static_chain,
        static_chain_res_count=query.static_chain_res_count,
        static_resn=query.static_resn,
        static_resi=query.static_resi,
        static_pos=query.static_pos,
        static_name=query.static_name,
        static_order=query.static_order,
        model_count=query.model_count,
        model=query.model,
        name=query.name,
        head=query.head,
        genes=query.genes,
        genus=query.genus,
        host=query.host,
        deposition=query.deposition,
        method=query.method,
        resolution=query.resolution,
        bio=query.bio,
        bio_list=query.bio_list,
        bio_count=query.bio_count,
        chain_count=query.chain_count,
        pdb_id=query.pdb_id,
        distance=query.distance,
        d1_X_C=query.d1_X_C,
        d2_X_CA=query.d2_X_CA,
        d3_X_N=query.d3_X_N,
        d4_CA_C=query.d4_CA_C,
        d5_CA_CA=query.d5_CA_CA,
        d6_CA_N=query.d6_CA_N,
        d7_N_C=query.d7_N_C,
        d8_N_CA=query.d8_N_CA,
        d9_N_N=query.d9_N_N,
        nterm=query.nterm,
        x_name=query.x_name,
        distances=query.distances,
        poly_ids=query.poly_ids,
        dynamic_chain_poly_id=query.dynamic_chain_poly_id,
        static_chain_poly_id=query.static_chain_poly_id,
        chain_poly_id_eq=query.chain_poly_id_eq,
        # inter=query.inter, #TODO remove comment for actual set
        dist_count=query.dist_count,
        inter_count=query.inter_count,
        inter=query.inter,
        hit_count=query.hit_count,
        cys_count=query.cys_count,
        job_id=str(entry),  # query.job_id, #TODO change this to be proper, good enough for now?
        timestamp=query.timestamp)
    return alchemy_sum


def submit_summary(summary, out_session):
    out_session.add(summary)
    try:
        out_session.commit()
        return True
        # print("\tCommitted entry")
    except exc.IntegrityError as e:
        out_session.rollback()
        print("\tsession.rollback() - {}".format(e))
        return False


# CONSTRUCTION ZONE - IN PROGRESS
def progressive_filter(entry, subquery):
    key = entry[0]
    term = entry[1]

    if key == 'distance':
        ssubquery = subquery.filter(DCSummary.distance > term)
    # TODO needs to throw error if undefined, or handle otherwise
    return subquery


def distance_series(query_file, db_input, db_output, cutoff, t1, t2, t3, dyn_name, poly_eq):
    db_in = create_engine(db_input)
    db_out = create_engine(db_output)
    DB_in_Session = sessionmaker(db_in)
    DB_out_Ssession = sessionmaker(db_out)
    in_session = DB_in_Session()
    out_session = DB_out_Ssession()

    # Query = in_session.query(DCSummary)
    start_time = time.time()
    parametered_query = in_session.query(Sub2DistanceUnique) #DIDN'T CAHNGE THIS LAST TIME
    print("Time (s): ", time.time() - start_time, "\n")

    distance_range = []
    parameter_range = []
    # TODO add option to SORT individual distances on DISTANCE table, THEN ONLY FIRST ENTRY, skip rest
    # parametered_query = Query.filter(DCSummary.abort == False)
    # parametered_query = parametered_query.filter(DCSummary.bio==1)

    parametered_count = parametered_query.count()
    print('Total count: ', parametered_count)
    print("Time (s): ", time.time() - start_time, "\n")

    print("Query parameters", t1, t2, t3, dyn_name)

    if dyn_name is False: #How does this behave
        pass
    else:
        parametered_query = parametered_query.filter(Sub2DistanceUnique.dynamic_name == dyn_name) #x_name for NZ, dynamic_name for Nt?
        #parametered_query = parametered_query.filter(CompDistanceUnique.dynamic_name == x_name) #x_name for NZ, dynamic_name for Nt?
        #corrected dynamic_name to x_name to handle N-terminus edge-case

    parametered_query = parametered_query.filter(Sub2DistanceUnique.min_type == t1)

    if poly_eq is False:
        pass
    else:
        parametered_query = parametered_query.filter(Sub2DistanceUnique.poly_id_eq == poly_eq)

    if t2 is False:
        pass
    else:
        parametered_query = parametered_query.filter(Sub2DistanceUnique.cnt_star == t2)
    if t3 is False:
        pass
    else:
        parametered_query = parametered_query.filter(Sub2DistanceUnique.cnt_pideq == t3)

    # parametered_query = parametered_query.filter(CompDistance.min_type == 'A') #all
    # parametered_query = parametered_query.filter(CompDistance.min_type == 'I0').filter(CompDistance.cnt_star == 2).filter(CompDistance.cnt_pideq==1) #intra
    # parametered_query = parametered_query.filter(CompDistance.min_type == 'I1E0')#.filter(CompDistance.cnt_pideq == 2) #hetero/inter
    # parametered_query = parametered_query.filter(CompDistance.min_type == 'I1E1') #homo/inter

    parametered_query = parametered_query.order_by(Sub2DistanceUnique.d7_N_C) #d1_X_C
    #distance to d1_X_C
    parametered_count = parametered_query.count()
    print('Total count: ', parametered_count)
    print("Time (s): ", time.time() - start_time, "\n")
    distance_range.append(("total", parametered_count))
    for i in range(0, cutoff):
        # distance_query = parametered_query.filter(DCSummary.distance < i/10)
        distance_query = parametered_query.filter(Sub2DistanceUnique.d7_N_C < i / 10) #d1_X_C instead of distance
        distance_count = 0
        distance_count = distance_query.count()
        distance_range.append((i, distance_count))
        print(i / 10, distance_count)
    print("Time (s): ", time.time() - start_time, "\n")

    """ #This crashes
    start_time = time.time()
    decrement_query = parametered_query
    for i in reversed(range(0,501)):
        decrement_query2 = decrement_query.filter(DCSummary.distance < i/10)
        del decrement_query
        decrement_query = decrement_query2
        distance_count = decrement_query.count()
        distance_range.append((i, distance_count))
        print('Cutoff: ', i / 10, 'Count: ', distance_count)
    print("Time (s): ", time.time() - start_time, "\n")
    """

    in_session.close()
    out_session.close()
    return distance_range


# END OF CONSTRUCTION ZONE

if __name__ == "__main__":

    """
    for i in range(0, 70):
        print('DELETE FROM tempdistance WHERE pdb_id like "%%";')
        print("INSERT INTO tempdistance SELECT *,_rowid_ from distance LIMIT 1000000 OFFSET {};".format(1000000*i))
        print("INSERT INTO subdistance SELECT *, min(distance),'A',_rowid_ from tempdistance GROUP BY pdb_id;")
        print("INSERT INTO subdistance SELECT *, min(distance),'I0', _rowid_ from tempdistance WHERE inter is 0 GROUP BY pdb_id;")
        print("INSERT INTO subdistance SELECT *, min(distance),'I1E0',_rowid_ from tempdistance WHERE inter is 1 and not poly_id_eq GROUP BY pdb_id;")
        print("INSERT INTO subdistance SELECT *, min(distance), 'I1E1',_rowid_ from tempdistance WHERE inter is 1 and poly_id_eq GROUP BY pdb_id;")
    """

    """
        print("INSERT INTO subdistance SELECT *, min(distance),'N',_rowid_ from tempdistance WHERE d7_N_C not NULL and nterm GROUP BY pdb_id;")
        print("INSERT INTO subdistance SELECT *, min(distance),'NI', _rowid_ from tempdistance WHERE inter is 0 and d7_N_C not NULL and nterm GROUP BY pdb_id;")
        print("INSERT INTO subdistance SELECT *, min(distance),'NIE',_rowid_ from tempdistance WHERE inter is 1 and d7_N_C not NULL and nterm and not poly_id_eq GROUP BY pdb_id;")
        print("INSERT INTO subdistance SELECT *, min(distance), 'NIE',_rowid_ from tempdistance WHERE inter is 1 and d7_N_C not NULL and nterm and poly_id_eq GROUP BY pdb_id;")
    """

    # main(query_file, db_input, db_output, log_file_csv)

    query_file = "C:\\Users\\Arne\\Desktop\\disCrawl\\Queries\\output\\query_tab.txt"
    #db_input = "sqlite:///C://Users//Arne//Desktop//Thesis//Chapter 3//disCrawl stats ch3//disCrawl_XtoC_redone.db"
    db_input = "sqlite:///C://Users//Arne//Desktop//Thesis//Chapter 3//disCrawl stats ch3//disCrawl-Nt-to-C.db"

    db_output = "sqlite:///C://Users//Arne//Desktop//disCrawl_filtered.db"
    # log_file = "C:\\Users\\Arne\\Desktop\\disCrawl\\Queries\\disCrawl_query_log.txt"
    log_file_csv = "C:\\Users\\Arne\\Desktop\\Thesis\\Chapter 3\\disCrawl stats ch3\\disCrawl_query_results.txt"

    """
    
    input_queue = [("I0",False,False), ("I0",1,False),("I0",2,False),("I0",2,1),("I0",2,2),("I0",3,False),
                   ("I1E0", False, False), ("I1E0", 1, False), ("I1E0", 2, False), ("I1E0", 2, 1), ("I1E0", 2, 2),
                   ("I1E0", 3, False),("I1E1", False, False), ("I1E1", 1, False), ("I1E1", 2, False), ("I1E1", 2, 1), ("I1E1", 2, 2),
                   ("I1E1", 3, False)]
    """
    """
    input_queue = [("NI", False, False, "NZ", False), ("NI", 1, False, "NZ", False), ("NI", 2, False, "NZ", False), ("NI", 2, 1, "NZ", False), ("NI", 2, 2, "NZ", False),
                   ("NI", 3, False, "NZ", False),
                   ("NIE", False, False,"NZ", 0), ("NIE", 1, False,"NZ", 0), ("NIE", 2, False,"NZ", 0), ("NIE", 2, 1,"NZ", 0), ("NIE", 2, 2,"NZ", 0),
                   ("NIE", 3, False,"NZ", 1), ("NIE", False, False,"NZ", 1), ("NIE", 1, False,"NZ", 1), ("NIE", 2, False,"NZ", 1), ("NIE", 2, 1,"NZ", 1),
                   ("NIE", 2, 2,"NZ", 1),
                   ("NIE", 3, False,"NZ", 1), ("NI", False, False, "OH", False), ("NI", 1, False, "OH", False), ("NI", 2, False, "OH", False), ("NI", 2, 1, "OH", False), ("NI", 2, 2, "OH", False),
                   ("NI", 3, False, "OH", False),
                   ("NIE", False, False,"OH", 0), ("NIE", 1, False,"OH", 0), ("NIE", 2, False,"OH", 0), ("NIE", 2, 1,"OH", 0), ("NIE", 2, 2,"OH", 0),
                   ("NIE", 3, False,"OH", 1), ("NIE", False, False,"OH", 1), ("NIE", 1, False,"OH", 1), ("NIE", 2, False,"OH", 1), ("NIE", 2, 1,"OH", 1),
                   ("NIE", 2, 2,"OH", 1),
                   ("NIE", 3, False,"OH", 1)]

    """
    """
    input_queue = [("N", False, False, "NZ", False), ("NI", False, False, "NZ", False), ("NIE", False, False, "NZ", 0),
                   ("NIE", False, False, "NZ", 1)
        , ("N", False, False, "OH", False), ("NI", False, False, "OH", False), ("NIE", False, False, "OH", 0),
                   ("NIE", False, False, "OH", 1)]
    """
    
    input_queue = [("N", False, False, False, False), ("NI", False, False, False, False), ("NIE", False, False, False, 0),
                   ("NIE", False, False, False, 1)
        , ("N", False, False, "N", False), ("NI", False, False, "N", False), ("NIE", False, False, "N", 0),
                   ("NIE", False, False, "N", 1)]

    """
    input_queue = [("NI", False, False, False, False), ("NI", 1, False, False, False), ("NI", 2, False, False, False), ("NI", 2, 1, False, False), ("NI", 2, 2, False, False),
                   ("NI", 3, False, False, False),
                   ("NIE", False, False,False, 0), ("NIE", 1, False,False, 0), ("NIE", 2, False,False, 0), ("NIE", 2, 1,False, 0), ("NIE", 2, 2,False, 0),
                   ("NIE", 3, False,False, 1), ("NIE", False, False,False, 1), ("NIE", 1, False,False, 1), ("NIE", 2, False,False, 1), ("NIE", 2, 1,False, 1),
                   ("NIE", 2, 2,False, 1),
                   ("NIE", 3, False,False, 1), ("NI", False, False, "N", False), ("NI", 1, False, "N", False), ("NI", 2, False, "N", False), ("NI", 2, 1, "N", False), ("NI", 2, 2, "N", False),
                   ("NI", 3, False, "N", False),
                   ("NIE", False, False,"N", 0), ("NIE", 1, False,"N", 0), ("NIE", 2, False,"N", 0), ("NIE", 2, 1,"N", 0), ("NIE", 2, 2,"N", 0),
                   ("NIE", 3, False,"N", 1), ("NIE", False, False,"N", 1), ("NIE", 1, False,"N", 1), ("NIE", 2, False,"N", 1), ("NIE", 2, 1,"N", 1),
                   ("NIE", 2, 2,"N", 1),
                   ("NIE", 3, False,"N", 1)]

    """
    cutoff = 501

    result_queue = []
    for el in input_queue:
        result_queue.append((el, distance_series(query_file, db_input, db_output, cutoff, *el)))
    print(result_queue)

    results_reordered = [[], []]

    results_reordered[1].append("total")
    for x in range(0, cutoff):
        results_reordered[1].append(x)
    for el in result_queue:
        results_reordered[0].append(el[0])
        results_reordered.append([])
        for x, y in el[1]:
            results_reordered[len(results_reordered) - 1].append(y)
    print(results_reordered)

    with open("C:\\Users\\Arne\\Desktop\\Thesis\\Chapter 3\\disCrawl stats ch3\\outfile-sub-20201230-corrected.txt", "w+") as outf:

        for i in range(0, cutoff + 1):
            line = ""
            if i == 0:
                line = line + "Cutoff"
                for el in results_reordered[0]:
                    line = line + "|" + str(el)
                print(line)
                outf.write(line + "\n")
                line = ""
            j = 0
            for el in results_reordered:
                if j == 0:
                    pass
                else:
                    line = line + "|" + str(el[i])
                j += 1
            print(line)
            outf.write(line.lstrip("|") + "\n")