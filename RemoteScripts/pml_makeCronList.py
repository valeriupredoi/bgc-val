#!/usr/bin/python


from ConfigParser import ConfigParser
import os

package_directory = os.path.dirname(os.path.abspath(__file__))
headerline =  "\n\n####"


def checkConfig(Config,debug=False):
	"""
	If it's a string, it opens it as a ConfigParser.
	"""
	if type(Config) == type('string'):
		if debug: print "Reading", Config
		Config1 = ConfigParser()
		Config1.read(Config)
		Config = Config1
	return Config

def parseList(Config,section,option):
	"""
	This tool loads an string from config file and returns it as a list.
	"""
	Config = checkConfig(Config)
	try:	list1 = Config.get(section, option)
	except: 
		print "No option ",option," in section: ",section
		return ''
		
	list1 = list1.replace(',', ' ')
	list1 = list1.replace('  ', ' ')
	list1 = list1.replace('\'', '')
	list1 = list1.replace('\"', '')	
	while list1.count('  ')>0: 
		list1 = list1.replace('  ', ' ')
	return list1.split(' ')




def addCompare(txt=''):
	txt += headerline
	txt += "\n# analysis_compare"
	txt += "\n1 8   * * * /users/modellers/ledm/workspace/ukesm-validation/RemoteScripts/pmlcron_compare.sh >  /users/modellers/ledm/ImagesFromJasmin/cron-compare.log 2>&1"
	return txt


def loadJobs(cp):
	cp = checkConfig(cp)
	jobs = cp.options('jobs')
	jobs = {j:True for j in jobs} #remove duplicates
	jobs = sorted(jobs.keys())
	jobs.reverse()
	return jobs


def hourcheck(hour):
	if 0 <= hour < 24:
		return hour
	while hour > 23:
		hour -= 24
		if hour <0: assert 0
	return hour

def addMassjobs(configfn='',startinghour = 18):
	cp = checkConfig(configfn,)
	jobs = loadJobs(cp)
	
	massscripts = ''
	massscripts += headerline	
	massscripts += "\n# DownloadFromMass"	
	for i, job in enumerate(jobs):
		hour = hourcheck(startinghour+i)

		options = 	parseList(cp,'jobs',job)
		if 'standard' in options:	
			massscripts	+= '\n1  '+str(int(hour))+\
					'   * * * /users/modellers/ledm/workspace/ukesm-validation/RemoteScripts/pmlcron_mass.sh '+\
					job + ' >  /users/modellers/ledm/ImagesFromJasmin/cron-mass_'+\
					job + '.log 2>&1'
		
		if 'MLD' in options:	
			massscripts	+= '\n30 '+str(int(hour))+\
					'   * * * /users/modellers/ledm/workspace/ukesm-validation/RemoteScripts/pmlcron_mass_MLD.sh '+\
					job + ' >  /users/modellers/ledm/ImagesFromJasmin/cron-mass_mld_'+\
					job + '.log 2>&1'
	return massscripts	

def addSci1jobs(configfn='',startinghour = 21):
	cp = checkConfig(configfn,)
	jobs = loadJobs(cp)
	
	scijobs = ''
	scijobs += headerline
	scijobs += "\n# Run evaluation"			
	for i, job in enumerate(jobs):
		hour = hourcheck(startinghour+i)

		options = 	parseList(cp,'jobs',job)
		if 'standard' in options:	
			scijobs	+= '\n1  '+str(int(hour))+\
					'   * * * /users/modellers/ledm/workspace/ukesm-validation/RemoteScripts/pmlcron_sci1.sh '+\
					job + ' >  /users/modellers/ledm/ImagesFromJasmin/cron-sci1_'+\
					job + '.log 2>&1'
		
		if 'physics' in options:	
			scijobs	+= '\n30 '+str(int(hour))+\
					'   * * * /users/modellers/ledm/workspace/ukesm-validation/RemoteScripts/pmlcron_sci1_physOnly.sh '+\
					job + ' >  /users/modellers/ledm/ImagesFromJasmin/cron_sci1_phys_'+\
					job + '.log 2>&1'
	return scijobs	

def runNow(configfn=''):
	#####
	# this script 

	cp = checkConfig(configfn,)
	jobs = loadJobs(cp)
	scijobs = ''
	scijobs += "\n# Run evaluation\n"			
	for i, job in enumerate(jobs):
		options = 	parseList(cp,'jobs',job)
		if 'standard' in options:	
			scijobs	+= '/users/modellers/ledm/workspace/ukesm-validation/RemoteScripts/pmlcron_mass.sh '+job + '; '
		
		if 'MLD' in options:	
			scijobs	+= '/users/modellers/ledm/workspace/ukesm-validation/RemoteScripts/pmlcron_mass_MLD.sh '+job + '; '
			
			
	for i, job in enumerate(jobs):
		options = 	parseList(cp,'jobs',job)
		if 'standard' in options:	
			scijobs	+= '/users/modellers/ledm/workspace/ukesm-validation/RemoteScripts/pmlcron_sci1.sh '+job + '; '
		
		if 'physics' in options:	
			scijobs	+= '/users/modellers/ledm/workspace/ukesm-validation/RemoteScripts/pmlcron_sci1_physOnly.sh '+job + '; '
		
	return scijobs	
	
	
		
def main():
	fn = str(package_directory+'/jobids_config.ini')

	crontxt = ''

	crontxt+= headerline+"\nMAILTO=\"\""
	
	crontxt+= headerline+"\n# This crontable was made by: "+str(os.path.abspath(__file__))
	
	crontxt+= headerline+"\n# using the file: "+str(os.path.abspath(fn))
	

	crontxt += addCompare()
	
	crontxt += addMassjobs(fn)
	
	crontxt += addSci1jobs(fn)
	
	crontxt+= headerline
	print crontxt
	save =True
	if save:
		cfn = 'crons/crontab.txt'
		f = open(cfn,'w')
		f.write(crontxt)
		f.close()
		print "saved as :", cfn
		print "install with: \ncrontab "+cfn
	
	print runNow(fn)
	
main()
