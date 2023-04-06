@echo off
setlocal EnableDelayedExpansion

set INPUT_FOLDER=C:\Users\worms\NGS Data\2022.06.07_drift_seq\90-666155004b\00_fastq\NNK\dada2_test\filtered
set OUTPUT_FOLDER="%INPUT_FOLDER%\merged_reads"
set INPUT_FOLDER="%INPUT_FOLDER%"

if not exist %OUTPUT_FOLDER% mkdir %OUTPUT_FOLDER%

for %%f in (%INPUT_FOLDER%\*_F_filt.fastq.gz) do (
    set READ1=%%f
	set "READ2=!READ1:_F_filt=_R_filt!"
    flash -o merged -M 110  -d !OUTPUT_FOLDER! "!READ1!" "!READ2!"
)


endlocal