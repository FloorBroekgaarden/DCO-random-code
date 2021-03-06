import cv2
import os
import string
import os
import moviepy.video.io.ImageSequenceClip
from PIL import Image
import glob


def makeMovie_DCOparam(DCOtype='BHNS', whichParam='Initial', fps=.4, duration=300):
	'''
	whichParam = 'Initial' or 'Final'
	fps=0.4, frames per second
	
	'''

	nModels=15
	BPSnameslist = list(string.ascii_uppercase)[0:nModels] 	

	image_folder = '/Users/floorbroekgaarden/Projects/BHNS_project/PlottingScripts/3_DCO-Population/'+DCOtype +'_DCOpropertiesModelVariations/'
	# images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
	images = []
	for ind_m, bps_model in enumerate(BPSnameslist):
	    images.append(image_folder +  whichParam + 'Param_'+ bps_model +'.png')


	image_files = images
	clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
	clip.write_videofile(image_folder+'movie_'+ whichParam + 'Param_' + DCOtype + '.mp4')


	# make also gif:
 
	# Create the frames
	frames = []
	# imgs = glob.glob("*.png")
	for i in images:
	    new_frame = Image.open(i)
	    frames.append(new_frame)
	 
	# Save into a GIF file that loops forever
	frames[0].save(image_folder+'gif_'+ whichParam + 'Param_' + DCOtype + '.gif', format='GIF',
	               append_images=frames[1:],
	               save_all=True,
	               duration=duration, loop=0)



	return 







def makeMovie_ObservedPDFsSFRDs(DCOtype='BHNS', bps_model='A', fps=.4, duration=300):
	'''
	whichParam = 'Initial' or 'Final'
	fps=0.4, frames per second
	
	'''




	GSMFs = ['Panter et al. (2004) Single', 'Furlong et al. (2015) Single', 'Furlong et al. (2015) Double']
	MZs   = [ 'Langer et al. (2006)'      , 'Langer et al. +offset (2006)', 'Ma et al. (2015)']
	SFRs  = ['Madau et al. (2014)'         ,'Strolger et al. (2004)',     'Madau et al. (2017)']


	MSSFRnameslist = []
	MSSFRnameslist.append('000') # add phenomenological 



	for ind_SFR, SFR in enumerate(SFRs):
		ind_x = ind_SFR + 1
		for ind_GSMF, GSMF in enumerate(GSMFs):
			ind_y = ind_GSMF + 1
			for ind_MZ, MZ in enumerate(MZs):
				ind_z = ind_MZ + 1
				
				MSSFRnameslist.append('%s%s%s'%(ind_x, ind_y, ind_z))



                
  
	image_folder = '/Users/floorbroekgaarden/Projects/BHNS_project/PlottingScripts/7_Discussion/DistributionsPlots'+DCOtype +'/'                 

	images = []
	for ind_m, SFRD_model in enumerate(MSSFRnameslist):
		images.append(image_folder +   'observed_' + bps_model + SFRD_model + '.png')


	image_files = images
	clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
	clip.write_videofile(image_folder+'movie_'+ 'ObservedDistributions_' + bps_model + '_' + DCOtype + '.mp4')


	# make also gif:
 
	# Create the frames
	frames = []
	# imgs = glob.glob("*.png")
	for i in images:
	    new_frame = Image.open(i)
	    frames.append(new_frame)
	 
	# Save into a GIF file that loops forever
	frames[0].save(image_folder+'gif_'+ 'ObservedDistributions_' + bps_model + '_' + DCOtype +  '.gif', format='GIF',
	               append_images=frames[1:],
	               save_all=True,
	               duration=duration, loop=0)



	return 



makeMovie_Observed=False
makeMovie_DCOproperties=True 



if makeMovie_Observed==True:
	makeMovie_ObservedPDFsSFRDs(DCOtype='BHNS', bps_model='A', fps=.4, duration=300)

# if makeMovie_DCOproperties==True:
# 	makeMovie_DCOparam(DCOtype='BHNS', whichParam='Initial', fps=.4, duration=300)
if makeMovie_DCOproperties==True:
	makeMovie_DCOparam(DCOtype='BHNS', whichParam='Final', fps=.4, duration=300)
