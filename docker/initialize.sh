#base="/Users/nanda/Documents/EclipseWorkspace/ATACpipeline"
#docker run -it  --volume $base"/pipeline-chromatin-accessibility:/pipeline-chromatin-accessibility" --volume $base"/data":/data docker_atac_platform

#!/bin/bash
docker build -t docker_atac_platform .
docker run -it  -d --volume pipeline-chromatin-accessibility:/pipeline-chromatin-accessibility --volume pipeline-chromatin-accessibility/data:/data docker_atac_platform
