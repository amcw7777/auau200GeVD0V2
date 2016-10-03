```bash
mkdir myAnalysis
cd myAnalysis

# Replace address below with your own fork if you have one
git clone git@github.com:amcw7777/auau200GeVD0V2.git
# Link all needed code under one StRoot directory:
mkdir StRoot
sh makeLinks.sh

# Compile
starver SL16d
cons
root4star -l -b -q -x runPicoD0AnaMaker.C\(\"$FILELIST\",\"OUTPUTFILENAME.root\"\)
```
