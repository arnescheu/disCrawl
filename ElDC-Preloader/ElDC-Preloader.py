from pymol import stored
from io import BytesIO
from io import StringIO
import os
import csv
import requests
import shutil

class dp: #dummy class for default parameters
	#location of this script for reinitialization
	script_location="ElDC-Preloader-Visualization.py"
	
	#csv input and hitlist parameters
	hitlist_path="dev\\PyMOL\\pdbIDs_hit.txt"
	csv_pdbID_path="\\dev\\PyMOL\\pdbIDs.txt"
	csv_delimiter="|"
	
	#structure file parameters
	obj_id_str='str'
	file_type_str='.cif'
	local_DB_str="pdb\\data\\structures\\divided\\mmCIF_as"
	local_DB_asym="pdb\\data\\structures\\divided\\mmCIF_ex"

	local_repo_str="dev\\PyMOL\\preloader\\structures\\cif\\"
	web_URL_str="https://files.rcsb.org/download/"
	#electron density file parameters
	obj_id_2mFo='2mFo'
	obj_id_mFo='mFo'
	file_type_2mFo=".ccp4"
	file_type_mFo="_diff.ccp4"
	local_DB_2mFo=None
	local_DB_mFo=None
	local_repo_2mFo="dev\\PyMOL\\preloader\\coordinates\\2mFo"
	local_repo_mFo="dev\\PyMOL\\preloader\\coordinates\\mFo"
	web_URL_map="http://www.ebi.ac.uk/pdbe/coordinates/files/"	
	#save to local repository after loading?
	save=True
	
	#isomesh parameters
	layers=3
	step=1
	start=0.5
	carve=2.5
	color_range=["blue","slate","lightblue","white"]

###LOADING SCRIPT###
def preload(pdbID_list, file_type=".cif", object_identifier="str", local_repository=dp.local_DB_str, local_DB="pdb\\data\\structures\\divided\\mmCIF_ex", web_URL="https://files.rcsb.org/download/", save=dp.save):
	"""Tries to load structures, from a database copy (in divided folder structure), or a flat local repository, or the web.
	"""	
	if isinstance(pdbID_list,str): pdbID_list=[pdbID_list]
	for pdbID in pdbID_list:	
		#loading object if new
		local_DB = dp.local_DB_str
		if ".asym" in pdbID:
			pdbID = pdbID.replace(".asym","")
			local_DB = dp.local_DB_asym
		if not (pdbID+"_"+object_identifier) in cmd.get_object_list('all'):
			if not local_repository is None:
				URL = os.path.join(local_repository,pdbID+file_type)
				if object_from_URL(URL,pdbID+"_"+object_identifier): 
					if file_type == ".cif":
						print_mmCIF_name(URL)
					continue
			if not local_DB is None:
				URL = os.path.join(local_DB,pdbID[1:3],pdbID+file_type)
				if object_from_URL(URL, pdbID+"_"+object_identifier): 
					if file_type == ".cif":
						print_mmCIF_name(URL)
					continue
			if not web_URL is None:
				#AS20180327: Changed to save file directly, not through pymol, to avoid degradation.
				#However, now it needs a local repository to save the URL to
				#save=False
				r=requests.get(web_URL+pdbID+file_type)
				if not r:
					print("Error "+str(r.status_code)+" in request from "+web_URL+pdbID+file_type)
				else:
					print("Loading from web: "+web_URL+pdbID+file_type)
					with open(os.path.join(local_repository,pdbID+file_type), "wb") as outfile:			
						shutil.copyfileobj(BytesIO(r.content), outfile)
					with open(os.path.join(local_repository,pdbID+file_type), "r") as outfile:			
						cmd.load(os.path.join(local_repository,pdbID+file_type),pdbID+"_"+object_identifier)
						#pdb_file.write(r.read())
					if not save:os.remove(os.path.join(local_repository,pdbID+file_type))
				continue
				#if object_from_URL(web_URL+pdbID+file_type,pdbID+"_"+object_identifier):
					#if save:
					#	cmd.save(os.path.join(local_repository,pdbID+file_type),pdbID+"_"+object_identifier)
					#	print("Saved to %s" %os.path.join(local_repository,pdbID+file_type))
					#continue
			print("Failed to load %s for %s" %(object_identifier,pdbID))
		else: print("%s_%s is preloaded" %(pdbID,object_identifier))
	
def object_from_URL(URL, object):
	"""Creates a PyMOL object (with the given name) from the given URL. Returns False if it fails"""
	print("Attempting to load %s from %s" %(object,URL))
	try:
		cmd.load(URL,object)
		print("%s loaded" %object)
		return True
	except CmdException:
		print("Failed to load %s from %s" %(object,URL))
		return False
		
def buffer_all():
	csv_set_pdbIDs()
	global pdbID_global_list
	
	local_repo_str=dp.local_repo_str
	for pdbID in pdbID_global_list:
		local_DB=dp.local_DB_str
		if ".asym" in pdbID:
			local_DB = dp.local_DB_asym
		try:
			shutil.copyfile(os.path.join(local_DB,pdbID[1:3],pdbID+".cif"), os.path.join(local_repo_str,pdbID+".cif"))
		except IOError as e:
			print(e)
	csv_set_pdbIDs()
	
def print_mmCIF_name(URL):
	with open(URL, "r") as mmCIF_file:
		for line in mmCIF_file:
			if "_struct.title" in line:
				print(line.replace("_struct.title","").strip())
				continue
			if "_struct.pdbx_descriptor" in line:
				print(line.replace("_struct.pdbx_descriptor","").strip())
				return True
		return False

###LOAD STRUCTURES WITH ISOMESH###
def load_map(pdbID_list,clear=True,carve_selector='all'):
	"""loads structure with mesh and applies colors"""
	if isinstance(pdbID_list,str): pdbID_list=[pdbID_list]
	if clear: cmd.delete("all")
	for pdbID in pdbID_list:
		preload([pdbID],file_type=dp.file_type_str, object_identifier=dp.obj_id_str, local_repository=dp.local_repo_str, local_DB=dp.local_DB_str, web_URL=dp.web_URL_str)
		preload([pdbID],file_type=dp.file_type_2mFo, object_identifier=dp.obj_id_2mFo, local_repository=dp.local_repo_2mFo, local_DB=dp.local_DB_2mFo, web_URL=dp.web_URL_map)
		preload([pdbID],file_type=dp.file_type_mFo, object_identifier=dp.obj_id_mFo, local_repository=dp.local_repo_mFo, local_DB=dp.local_DB_mFo, web_URL=dp.web_URL_map)
		isomash([pdbID], object_identifier=dp.obj_id_2mFo, carve_selector=carve_selector, layers=dp.layers, step=dp.step, start=dp.start, carve=dp.carve, color_range=dp.color_range)
		#except:print("Failed to isomash "+pdbID+"_2mFo")
		isomash([pdbID], object_identifier=dp.obj_id_mFo, carve_selector=carve_selector,layers=2, step=6, start=-3, carve=dp.carve, color_range=["green","red"])
		#except:print("Failed to isomash "+pdbID+"_mFo")
		isomashbow()

def isomash(pdbID_list=None, object_identifier='2mFo', carve_selector='all', layers=dp.layers, step=dp.step, start=dp.start, carve=dp.carve, color_range=dp.color_range):
	"""Applies isomesh layers."""
	print("Visualizing map for %s with %s start, %s carve, %s layers, %s step, color_range %s" %(pdbID_list, start, carve, layers, step, color_range))
	if pdbID_list is None:
		pdbID_list=[]
		for object in cmd.get_names():
			if "_str" in object:
				pdbID_list.append(object.split("_")[0])
	if layers>len(color_range): 
		layers=len(color_range)
		print("Limiting layers to color range. Reduce layers or add more colors.")
	if isinstance(pdbID_list,str): pdbID_list=[pdbID_list]
	print(pdbID_list[0]+" and "+carve_selector)
	for pdbID in pdbID_list:
		contour=start
		for level in range(0,layers):
			#cmd.isomesh("%s_%s_eld_%s"%(pdbID,object_identifier,level),"%s_%s"%(pdbID,object_identifier),contour,pdbID,carve=carve)
			try:
				cmd.isomesh("%s_%s_eld_%s"%(pdbID,object_identifier,level),"%s_%s"%(pdbID,object_identifier),contour,"%s_str and %s extend 5" %(pdbID,carve_selector),carve=carve)
			except:
				print("Failed to isomesh "+"%s_%s"%(pdbID,object_identifier)+" on "+"%s_str and %s" %(pdbID,carve_selector))
				break
			contour=contour+step
			cmd.color(color_range[level],"%s_%s_eld_%s"%(pdbID,object_identifier,level))
			
def isomashbow():
	"""Orients structure and applies colors"""
	cmd.hide("all")
	for object in cmd.get_names():
		if "_str" in object:
			cmd.show("sticks",object)
			util.cbc(selection=object)
			cmd.orient(selection=object)
			for ch in cmd.get_chains(object):
				cmd.color("red","br. last resn * and !hetatm and %s and chain %s" %(object,ch))
				cmd.color("blue","br. first resn * and !hetatm and %s and chain %s" %(object,ch))
			util.cnc(selection=object)
			cmd.orient(object)
		if "_eld" in object:
			cmd.show("mesh",object)
			
def next_map():
	"""Loads the next structure in pdbID_global_list with density map"""
	global pdbID_global_list
	try:
		load_map([pdbID_global_list[0]],clear=True)
		pdbID_global_list=pdbID_global_list[1:]
	except IndexError:
		print("End of list reached. Set list with set_pdbIDs(list) or csv_set_pdbIDs(filepath)")
				
#DISCRAWL VISUALIZATION
def visualize_disCrawl(cutoff=15,targetA_C=False,targetA_NZ=True):
	#"disCrawl": Lysines, Termini, Distances, Orientation
	for object in cmd.get_names():
		if "dis" in object:print("Ignore 'dis' object")
		else:
			cmd.hide("dots",object+" & resn hoh")
			cmd.show("cartoon",object)
			util.cbc(selection=object)
			cmd.color("grey",object+" & resn lys & n. C*")
			cmd.color("yellow",object+" & resn cys & n. C*")
			cmd.hide("lines",object)
			cmd.show("lines",object+" & resn lys")
			util.cnc(selection=object)
			cmd.orient(selection=object)
			for ch in cmd.get_chains(object):
				cmd.color("red","br. last resn * and !hetatm and /"+object+"//"+ch)
				cmd.color("blue","br. first resn * and !hetatm and /"+object+"//"+ch)
				cmd.show("lines","br. last resn * and !hetatm and /"+object+"//"+ch)
				cmd.show("lines","br. first resn * and !hetatm and /"+object+"//"+ch)
				cmd.select("CT","last n. C and !hetatm and /"+object+"//"+ch)
				cmd.select("prox_Lys","(resn lys  and /"+object+") in (bb. w. "+str(cutoff)+" of br. CT)")
				if targetA_C: cmd.select("prox_C","n. C & br. prox_Lys")
				if targetA_NZ: cmd.select("prox_NZ","n. NZ & br. prox_Lys")
				if targetA_NZ: cmd.distance('dis_NZ_'+object,'CT','prox_NZ')
				if targetA_C: cmd.distance('dis_C_'+object,'CT','prox_C')
				cmd.delete("CT")
				cmd.delete("prox_Lys")
				cmd.delete("prox_NZ")		
		
def load_dc(pdbID_list,clear=True): #load_disCrawl
	if isinstance(pdbID_list,str): pdbID_list=[pdbID_list]
	if clear: cmd.delete("all")
	for pdbID in pdbID_list:
		preload([pdbID],file_type=dp.file_type_str, object_identifier=dp.obj_id_str, local_repository=dp.local_repo_str, local_DB=dp.local_DB_str, web_URL=dp.web_URL_str)
	visualize_disCrawl(cutoff=15,targetA_C=False,targetA_NZ=True)
		
def next_dc():
	"""Loads the next structure in pdbID_global_list with density map"""
	global pdbID_global_list
	try:
		load_dc([pdbID_global_list[0]],clear=True)
		pdbID_global_list=pdbID_global_list[1:]
	except IndexError:
		print("End of list reached. Set list with set_pdbIDs(list) or csv_set_pdbIDs(filepath)")
		
###ROTATING THROUGH LIST OF IDs###
def set_pdbIDs(pdbID_list):
	"""Sets the current list of pdbIDs to rotate through from a list"""
	if isinstance(pdbID_list,str): pdbID_list=[pdbID_list]
	global pdbID_global_list
	pdbID_global_list=pdbID_list
	
def csv_set_pdbIDs(path=dp.csv_pdbID_path,delimiter=dp.csv_delimiter): #This loads pdbs seperated by line or delimiter.
	"""Sets the current list of pdbIDs to rotate through from the specified file"""
	with open(path, 'r') as dR:
		pdbIDs=[]
		input = list(csv.reader(dR, delimiter=delimiter))
		for line in input:
			pdbIDs.append(line[0])
			#for chain in line:
				#if chain == 0: pdbIDs.append(chain)
	global pdbID_global_list
	pdbID_global_list=pdbIDs
	
def next_str():
	"""Loads the next structure in pdbID_global_list without density map"""
	global pdbID_global_list
	try:
		cmd.delete("all")
		preload([pdbID_global_list[0]],dp.file_type_str,dp.obj_id_str,dp.local_repo_str,dp.local_DB_str,dp.web_URL_str)
		pdbID_global_list=pdbID_global_list[1:]
	except IndexError:
		print("End of list not loaded. Set list with set_pdbIDs(list) or csv_set_pdbIDs(filepath)")

###AUXILIARY###
###SAVING INTERESTING IDS###
def save_hit():
	"""Saves hit into a hitlist"""
	global hitlist
	try:
		for object in cmd.get_names():
			if not "dis" in object: 
				print("Writing %s to hitlist." %object)
				hitlist.write("%s\n" %object)
	except NameError:
		print("Initializing hitlist. Write to file with write_hits(path,mode).")
		hitlist = StringIO.StringIO()
		for object in cmd.get_names():
			if not "dis" in object: 
				print("Writing %s to hitlist." %object)
				hitlist.write("%s\n" %object)
	write_hits(path=(dp.hitlist_path).replace(".txt","-temp.txt"), mode="w+")
	
def write_hits(path=dp.hitlist_path, mode="w+"):
	"""Appends hitlist to file. Mode w: overwrite; mode a: append"""
	global hitlist
	try:
		with open(path, mode) as dR:
			dR.write(hitlist.getvalue()) 
	except NameError:
		print("Hitlist not initialized. Save with save_hit() or clear with clear_hitlist()")
			
def clear_hitlist():
	"""Clears hitlist"""
	global hitlist
	hitlist=StringIO.StringIO()
	print("Hitlist cleared")

def toggle_eld():
	"""Toggles electron density objects (by deleting them... can be slow)"""
	for object in cmd.get_names():
		if "_eld" in object: 
			cmd.delete(object)
		else: isomesh_visualization([object[0:4]],dp.layers,dp.step,dp.start,dp.carve,dp.color_range)
	
###EXPERIMENTAL FEATURES###
def csv_set_pdbIDs_sele(path=dp.csv_pdbID_path,delimiter=dp.csv_delimiter): #in construction
	"""Sets the current list of pdbIDs to rotate through from the specified file and includes additional information. Takes default pymol selection statements, tube seperated"""
	with open(path, 'r') as dR:
		pdbID_list=[]
		selection_list=[]
		reader = list(csv.DictReader(dR, delimiter=delimiter))
		for row in reader:
			pdbID=row['pdbID']
			selection=''
			for key in row:
				allowed=['chain','resi','resn','name','selection']
				if key in allowed: 
					if row[key] is None: print("ignore None at %s"%key)
					else: 
						if len(selection)==0:
							selection="%s %s" %(key, row[key])
						else: selection="%s and %s %s" %(selection, key, row[key])
				else: print("ignore %s. Allowed keys: %s" %(key,allowed))
				print(selection)
			pdbID_list.append(pdbID)
			selection_list.append(selection)

	global pdbID_global_list
	global sele_global_list

	pdbID_global_list=pdbID_list
	sele_global_list=selection_list
	
def next_map_sele():
	"""Loads the next structure in pdbID_global_list with density map"""
	global pdbID_global_list
	global sele_global_list
	try:
		load_with_isomesh([pdbID_global_list[0]],carve_selector=sele_global_list[0],clear=True)
		cmd.select("%s_hit"%pdbID_global_list[0],"%s_str and %s" %(pdbID_global_list[0],sele_global_list[0]))
		cmd.orient("%s_str and %s" %(pdbID_global_list[0],sele_global_list[0]))
		pdbID_global_list=pdbID_global_list[1:]
		d=len(sele_global_list)-len(pdbID_global_list)
		sele_global_list=sele_global_list[d:]
	except IndexError:
		print("End of list reached. Set list with set_pdbIDs(list) or csv_set_pdbIDs(filepath)")

def reinitialize_script(script=dp.script_location):
	"""reloads this script."""
	print("What you are doing is dangerous...\n Reinitializing...")
	print("Running "+script)
	cmd.run(script)
	
#HELLO WORLD#
print("Loaded structure and electron density loader from Arne HA Scheu")