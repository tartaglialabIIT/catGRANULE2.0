import numpy as np
import os
import pandas as pd
import subprocess
from multiprocessing import Pool


def list_jpg_files(directory):
    jpg_files = []
    for file in os.listdir(directory):
        if file.endswith(".jpg"):
            jpg_files.append(file)
    return jpg_files

def download_image(url):
    wget_cmd = ["wget", url, "-P", input_folder]
    result = subprocess.run(wget_cmd, stdout=subprocess.PIPE)
    return result

input_folder="./HPA/HPA_images/"
if os.path.isdir(input_folder)==False:
    os.mkdir(input_folder)

if os.path.isdir("./FINAL_RESULTS_HPA/")==False:
        os.mkdir("./FINAL_RESULTS_HPA/")
    
# Read the Supplementary table S3 from DeepPhase paper
table_S3=pd.read_excel("tableS3.xlsx",sheet_name="tableS3")
table_S3["image_name"]=table_S3.Image_url.str.split('/').str[-1]


num_processes = 24
    
    # Create a pool of workers
with Pool(num_processes) as pool:
    # Map the download function to the list of URLs
    results = pool.map(download_image, list(table_S3.Image_url))

for myurl in list(table_S3.Image_url):
    wget_cmd = ["wget", myurl, "-P", input_folder]
    result = subprocess.run(wget_cmd, stdout=subprocess.PIPE)

N=200

table_S3=table_S3.loc[:,["Image_url","image_name"]]
table_S3=table_S3.drop_duplicates()

print(len(table_S3))

num_chunks = len(list(set(table_S3.Image_url))) // N

# Loop over images sets
for i in range(num_chunks+1):
    start_idx_top = i * N
    end_idx_top = (i + 1) * N
    
    print(i,start_idx_top,end_idx_top)
    
    top_chunk = table_S3.iloc[start_idx_top:end_idx_top]
    # create a temporary folder to put top 200 images    
    tmp_dir="./tmp_input_final/"
    if os.path.isdir(tmp_dir)==False:
        os.mkdir(tmp_dir)
    
    # copy the images to the folders	
    for im in list(top_chunk.image_name):
        cp_cmd=["cp", input_folder+im, tmp_dir]
        result = subprocess.run(cp_cmd, stdout=subprocess.PIPE)
    
    # Run CellProfiler4
    cprof_cmd = ["cellprofiler", "-c", "-r", "-p", "cp4_pipeline.cppipe", "-i", tmp_dir, "-o", "FINAL_RESULTS_HPA/"]
    result = subprocess.run(cprof_cmd, stdout=subprocess.PIPE)

    # Rename the output files of CellProfiler
    out1_old="./FINAL_RESULTS_HPA/CP4Image.csv"
    out2_old="./FINAL_RESULTS_HPA/CP4cell.csv"
    out3_old="./FINAL_RESULTS_HPA/CP4RelateObjectsG.csv"
    out4_old="./FINAL_RESULTS_HPA/CP4RelateCellNucleus.csv"
    out5_old="./FINAL_RESULTS_HPA/CP4nucleus.csv"

    out1_new="./FINAL_RESULTS_HPA/CP4Image_"+str(i)+".csv"
    out2_new="./FINAL_RESULTS_HPA/CP4cell_"+str(i)+".csv"
    out3_new="./FINAL_RESULTS_HPA/CP4RelateObjectsG_"+str(i)+".csv"
    out4_new="./FINAL_RESULTS_HPA/CP4RelateCellNucleus_"+str(i)+".csv"
    out5_new="./FINAL_RESULTS_HPA/CP4nucleus_"+str(i)+".csv"

    
    ren_1=["mv" , out1_old, out1_new]
    result = subprocess.run(ren_1, stdout=subprocess.PIPE)
    ren_2=["mv" , out2_old, out2_new]
    result = subprocess.run(ren_2, stdout=subprocess.PIPE)
    ren_3=["mv" , out3_old, out3_new]
    result = subprocess.run(ren_3, stdout=subprocess.PIPE)
    ren_4=["mv" , out4_old, out4_new]
    result = subprocess.run(ren_4, stdout=subprocess.PIPE)
    ren_5=["mv" , out5_old, out5_new]
    result = subprocess.run(ren_5, stdout=subprocess.PIPE)

    # Remove the input folder
    rm_cmd=["rm", "-r",tmp_dir]
    result = subprocess.run(rm_cmd, stdout=subprocess.PIPE)
