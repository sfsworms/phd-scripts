@echo off
setlocal EnableDelayedExpansion

rem This calls the FLASH software...

set "INPUT_FOLDER=C:\Users\Sebastian Worms\ngs_data\2023 NGS Ale\90-933598625\sample_100k\filtered"
set "OUTPUT_FOLDER=%INPUT_FOLDER%\merged_reads"

echo Input folder: %INPUT_FOLDER%
echo Output folder: %OUTPUT_FOLDER%

if not exist "%OUTPUT_FOLDER%" mkdir "%OUTPUT_FOLDER%"

for %%f in ("%INPUT_FOLDER%\*F_filt.fastq.gz") do (
    set "READ1=%%f"
    set "READ2=!READ1:F_filt.fastq.gz=R_filt.fastq.gz!"
    
    echo Processing: !READ1!
    echo With: !READ2!
    flash -o merged -M 93 -d "!OUTPUT_FOLDER!" "!READ1!" "!READ2!" -z
)

endlocal

flash -o merged -M 93 -d "C:\Users\Sebastian Worms\ngs_data\2023 NGS Ale\90-933598625\00_fastq\merged_reads" "C:\Users\Sebastian Worms\ngs_data\2023 NGS Ale\90-933598625\00_fastq\filtered\pinholin-library_F_filt.fastq.gz" "C:\Users\Sebastian Worms\ngs_data\2023 NGS Ale\90-933598625\00_fastq\filtered\pinholin-library_R_filt.fastq.gz" -z