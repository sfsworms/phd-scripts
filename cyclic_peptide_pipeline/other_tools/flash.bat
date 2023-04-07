@echo off

setlocal EnableDelayedExpansion
setlocal EnableDelayedExpansion

rem This calls the FLASH software. INPUT_FOLDER should be a path to a folder with one fwr reads and the corresponding rv reads. The software assume forward reads ends in _F_filt.fastq.gz and the reverse reads in _R_filt.fastq.gz but have otherwise identical names. The naming scheme comes from dada2


set INPUT_FOLDER=C:\Users\worms\NGS Data\2022.06.07_drift_seq\90-666155004b\00_fastq\NNK\filtered\gen_1
set OUTPUT_FOLDER="%INPUT_FOLDER%\merged_reads"
set INPUT_FOLDER="%INPUT_FOLDER%"

if not exist %OUTPUT_FOLDER% mkdir %OUTPUT_FOLDER%

rem - M is the max overlap size. -o specify the output name, -d the output directory, -z gzip the output

for %%f in (%INPUT_FOLDER%\*_F_filt.fastq.gz) do (
    set READ1=%%f
	set "READ2=!READ1:_F_filt=_R_filt!"
    flash -o merged -M 93 -d !OUTPUT_FOLDER! "!READ1!" "!READ2!" -z
)

endlocal