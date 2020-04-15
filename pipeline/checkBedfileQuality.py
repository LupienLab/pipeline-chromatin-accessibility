import re
import shutil

#####################################################################
def checkBedFile(bed):
    with open(bed,'r') as BEDFH:
        content = BEDFH.readlines()
        firstLine = content[0]
        if(re.match(r"^[a-zA-Z]+[\s\t]+[a-zA-Z]+[\s\t]+[a-zA-Z]+", firstLine)):
            with open("tmpBed",'w') as TMPBEDFH:
                TMPBEDFH.write("".join(content[1:]))
                shutil.move("tmpBed", bed)
            
#####################################################################

