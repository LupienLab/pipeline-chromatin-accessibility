#docker build -t docker_atac_platform .
base="/Users/ankita/Documents/eclipse_workspace/ATAC_platform"
mkdir $base"/data"
cp $base"/pipeline-chromatin-accessibility/pipeline/Snakefile" $base"/data"
docker run -it  --volume $base"/pipeline-chromatin-accessibility:/pipeline-chromatin-accessibility" --volume $base"/data":/data docker_atac_platform
