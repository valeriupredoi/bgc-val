#!/usr/bin/ipython 

#
# Copyright 2015, Plymouth Marine Laboratory
#
# This file is part of the bgc-val library.
#
# bgc-val is free software: you can redistribute it and/or modify it
# under the terms of the Revised Berkeley Software Distribution (BSD) 3-clause license. 

# bgc-val is distributed in the hope that it will be useful, but
# without any warranty; without even the implied warranty of merchantability
# or fitness for a particular purpose. See the revised BSD license for more details.
# You should have received a copy of the revised BSD license along with bgc-val.
# If not, see <http://opensource.org/licenses/BSD-3-Clause>.
#
# Address:
# Plymouth Marine Laboratory
# Prospect Place, The Hoe
# Plymouth, PL1 3DH, UK
#
# Email:
# ledm@pml.ac.uk

import os



def AddtoFile(filepath,linenumber,text):
	"""	Adds "text" at line "linenumber" in the file "filepath"
	"""
	f = open(filepath, "r")
	contents = f.readlines()
	f.close()

	contents.insert(linenumber, text)

	f = open(filepath, "w")
	contents = "".join(contents)
	f.write(contents)
	f.close()

def locateLineNumber(filepath,linekey):
	"""	Locates the line above the first instance of "linekey" in the file "filepath."
	"""
	
	f = open(filepath, "r")
	contents = f.readlines()
	f.close()
	for l,line in enumerate(contents):
		if line.find(linekey)>=0:
			return l

def writeSideBar(filepath, href, text):
	"""	Addes a sidebar to the file "filepath."
	"""	
								
	newline = '\n\t\t\t\t\t\t\t\t<li><a href="#'+href+'" id="'+href+'-link" class="skel-layers-ignoreHref"><span class="icon fa-th">'+text+'</span></a></\n'
	linenumber = locateLineNumber(filepath, 'AddSideBar') -1
	print "Adding at ", linenumber,"line:\n",newline
	AddtoFile(filepath,linenumber,newline)



def writeDescription(filepath, text):
	"""	Addes a sidebar to the file "filepath."
	"""	
	newline = '\t\t\t\t\t\t\t<p>'+text+'</p>'
	linenumber = locateLineNumber(filepath, 'HeaderDescription') 
	print "Adding at ", linenumber,"line:\n",newline
	AddtoFile(filepath,linenumber,newline)


def fnToTitle(fn):
	return os.path.basename(fn).replace('.png', '').replace('_',' ')

def addImagesText(imagePath,title = ''):
	""" Creates the appropriate text to insert an image onto the page.
	"""
	f = open("html5/html5Assets/figure-template.html", "r")
	contents = f.readlines()
	f.close()
	if title == '': title = fnToTitle(imagePath)
	for l,line in enumerate(contents):
		if line.find('imagefilename')>=0:	contents[l] = contents[l].replace('imagefilename',imagePath)
		if line.find('PlotTitle')>=0:		contents[l] = contents[l].replace('PlotTitle',title)

	outtxt = '\n'.join(contents)		
	return outtxt

def AddSection(filepath,href,Title, Description='',Files=[]):
	"""	Addes a section to the file "filepath."
	"""
	#####
	# Add a link to the side bar
	writeSideBar(filepath, href, Title)

	##### 
	# Copy the template and add the images.
	f = open("html5/html5Assets/section-template.html", "r")
	contents = f.readlines()
	f.close()
	if type(Files) == type(['a','list',]):
		imagesTxt = '\n'.join([addImagesText(f) for f in Files])
		
	if type(Files) == type({'a':'dict',}):
		imagesTxt = '\n'.join([addImagesText(f,title=Files[f]) for f in sorted(Files.keys())])
			
	#####
	# Add Title, description and figures to the copied template
	for l,line in enumerate(contents):
		if line.find('href')>=0:	contents[l] = contents[l].replace('href',href)
		if line.find('Title')>=0:	contents[l] = contents[l].replace('Title',Title)
		if line.find('Description')>=0:	contents[l] = contents[l].replace('Description',Description)				
		if line.find('Figures')>=0:	contents[l+1] += imagesTxt
	
	#####
	# Convert the list into a string				
	outtxt = '\n'.join(contents)
	
	#####
	# Add this into the template file.
	linenumber = locateLineNumber(filepath, 'AddSectionHere') -1
	
	AddtoFile(filepath,linenumber,outtxt)	

#def AddImage(filepath, href, text):
#	"""	Addes a sidebar to the file "filepath."
#	"""	
#	newline = '<li><a href="#'+href+'" id="'+href+'-link" class="skel-layers-ignoreHref"><span class="icon fa-th">'+text+'</span></a></'
#	linenumber = locateLineNumber(filepath, 'AddSideBar')
#	print "Adding at ", linenumber,"line:\n",newline
#	AddtoFile(filepath,linenumber,newline)
	
	
