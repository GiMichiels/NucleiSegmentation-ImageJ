import sys
import os
from cytomine import Cytomine
from cytomine.models import *
from subprocess import call
from skimage import io
import numpy as np
from sldc import locator
from sldc.locator import affine_transform
from shapely.ops import cascaded_union
from argparse import ArgumentParser

def _upload_annotation(cytomine, img_inst, polygon, label=None, proba=1.0):
    image_id = img_inst.id
    # Transform in cartesian coordinates
    polygon = affine_transform(xx_coef=1, xy_coef=0, yx_coef=0, yy_coef=-1, delta_y=img_inst.height)(polygon)

    annotation = cytomine.add_annotation(polygon.wkt, image_id)
    if label is not None and annotation is not None:
        annotation_term = AnnotationTerm()
        annotation_term.annotation = annotation.id
        annotation_term.annotationIdent = annotation.id
        annotation_term.userannotation = annotation.id
        annotation_term.term = label
        annotation_term.expectedTerm = label
        annotation_term.rate = proba
        cytomine.save(annotation_term)
    return annotation

def download_images(images, destDir):
	# download the images
	for image in images:
		# url format: CYTOMINEURL/api/imageinstance/$idOfMyImageInstance/download
		url = cytomine_host+"/api/imageinstance/" + str(image.id) + "/download"
		filename = str(image.id) + ".tif"
		conn.fetch_url_into_file(url, destDir+"/"+filename, True, True) 

def remove_label(filename, label):
	result = filename.replace(label, '')
	return result

baseOutputFolder = "/dockershare/";

parser = ArgumentParser(prog="IJSegmentClusteredNuclei.py", description="ImageJ workflow to segment clustered nuclei")
parser.add_argument('--cytomine_host', dest="cytomine_host", default='http://localhost-core')
parser.add_argument('--cytomine_public_key', dest="cytomine_public_key", default="77af7d84-b737-4489-8864-d5ad93f4700b")
parser.add_argument('--cytomine_private_key', dest="cytomine_private_key", default="3ef1f34e-2a9e-4ff8-96ec-df8e57c7dfcd")
parser.add_argument("--cytomine_id_project", dest="cytomine_id_project", default="5378")
parser.add_argument("--ij_radius,", dest="radius", default="5")
parser.add_argument('--ij_threshold', dest="threshold", default="-0.5")
arguments, others = parser.parse_known_args(sys.argv)
radius = arguments.radius
threshold = arguments.threshold

#Cytomine connection parameters
cytomine_host=arguments.cytomine_host

id_project=arguments.cytomine_id_project

conn = Cytomine(arguments.cytomine_host, arguments.cytomine_public_key, arguments.cytomine_private_key, base_path = '/api/', working_path = '/tmp/', verbose= True)

current_user = conn.get_current_user()
user_job = current_user

job = conn.get_job(user_job.job)
# job=666

job = conn.update_job_status(job, status = job.RUNNING, progress = 0, status_comment = "Loading images...")

# Get the list of images in the project
image_instances = ImageInstanceCollection()
image_instances.project  =  id_project
image_instances  =  conn.fetch(image_instances)
images = image_instances.data()

# split the list of all images into a list of input images and 
# a list of mask images representing the ground truth data.

inputImages = []
masks = []
fileIDs = {}
for image in images:
	if "_lbl." in image.filename:
		masks.append(image)
	else:
		inputImages.append(image)
		head, tail = os.path.split(image.filename)
		fileIDs[tail] = image.id

# create the folder structure for the folders shared with docker 
jobFolder = baseOutputFolder + str(job.id) + "/"
#jobFolder = baseOutputFolder + str(job) + "/"
inDir = jobFolder + "in"
outDir = jobFolder + "out"
groundTruthDir = jobFolder + "ground_truth"

if not os.path.exists(inDir):
    os.makedirs(inDir)

if not os.path.exists(outDir):
    os.makedirs(outDir)

if not os.path.exists(groundTruthDir):
    os.makedirs(groundTruthDir)

# download the images
download_images(inputImages, inDir)

# call the image analysis workflow in the docker image
shArgs = "data/in data/out "+radius+" "+threshold + ""
job = conn.update_job_status(job, status = job.RUNNING, progress = 25, status_comment = "Launching workflow...")

command = "/fiji/run.sh " + shArgs

#command = "docker run --rm -v "+jobFolder+":/fiji/data neubiaswg5/nucleisegmentation-imagej " + shArgs
call(command,shell=True)	# waits for the subprocess to return

job = conn.update_job_status(job, status = job.RUNNING, progress = 40, status_comment = "Extracting polygons...")

for image in inputImages:
	file = str(image.id) + ".tif"
	path = outDir + "/" + file
	imageData = io.imread(path)

  	indexes = np.unique(imageData)

    	# locate polygons
    	locobj = locator.BinaryLocator()
    	objects = dict()

    	for i, index in enumerate(indexes):
		if index == 0:
		    continue
		mask = (imageData == index).astype(np.uint8) * 255
		polygons = [polygon[0].buffer(2.0) for polygon in locobj.locate(mask)]
		polygon = cascaded_union(polygons).buffer(-2.0)
		if not polygon.is_empty and polygon.area > 0:
		    objects[index] = polygon
	 	
	print("Found {} polygons in this images.".format(len(objects)))	
	slide = conn.get_image_instance(image.id)
	
	# upload
	for index, polygon in objects.items():
	        annotation = _upload_annotation(conn, slide, polygon, label=None)
	        if annotation:
			conn.add_annotation_property(annotation.id, "index", str(index))
	
job = conn.update_job_status(job, status = job.TERMINATED, progress = 50, status_comment =  "Cleaning up..")

# cleanup - remove the downloaded images 

for image in inputImages:
	file = str(image.id) + ".tif"
	path = inDir + "/" + file
	os.remove(path);

job = conn.update_job_status(job, status = job.TERMINATED, progress = 60, status_comment =  "Downloading ground truth data..")

# download ground truth data
for image in masks:
	url = arguments.cytomine_host+"/api/imageinstance/" + str(image.id) + "/download"
	inputImageFilename = remove_label(image.filename, "_lbl")
	head, tail = os.path.split(inputImageFilename)
	filename = str(fileIDs[tail]) + ".tif"
	conn.fetch_url_into_file(url, groundTruthDir+"/"+filename, True, True)

job = conn.update_job_status(job, status = job.TERMINATED, progress = 100, status_comment =  "Finished Job..")
