#!/bin/bash

# Define input and output directories
input_dir="/vol/storage/dest_linvilla_bams"
output_dir="/vol/storage/dest_linvilla_bams/dest_linvilla_bams_fixed"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through BAM files and add unique read groups
for bamfile in "$input_dir"/*.bam; do
    filename=$(basename "$bamfile")
    
    # Assign unique sample names based on filename
    case "$filename" in
        "US_Pen_Lin_1_2009-07-15.mel.bam") sample="Sample1" ;;
		"US_Pen_Lin_1_2009-11-16.mel.bam") sample="Sample2" ;;
		"US_Pen_Lin_1_2010-07-15.mel.bam") sample="Sample3" ;;
		"US_Pen_Lin_1_2010-11-16.mel.bam") sample="Sample4" ;;
		"US_Pen_Lin_1_2011-07-15.mel.bam") sample="Sample5" ;;
		"US_Pen_Lin_1_2011-10-15.mel.bam") sample="Sample6" ;;
		"US_Pen_Lin_1_2011-11-16.mel.bam") sample="Sample7" ;;
		"US_Pen_Lin_1_2012-07-16.mel.bam") sample="Sample8" ;;
		"US_Pen_Lin_1_2012-09-13.mel.bam") sample="Sample9" ;;
		"US_Pen_Lin_1_2012-10-10.mel.bam") sample="Sample10" ;;
		"US_Pen_Lin_1_2013-06-15.mel.bam") sample="Sample11" ;;
		"US_Pen_Lin_1_2014-06-27.mel.bam") sample="Sample12" ;;
		"US_Pen_Lin_1_2014-10-17.mel.bam") sample="Sample13" ;;
		"US_Pen_Lin_1_2015-07-15.mel.bam") sample="Sample14" ;;
		"US_Pen_Lin_1_2015-10-15.mel.bam") sample="Sample15" ;;
        *) echo "Skipping unknown BAM file: $filename"; continue ;;
    esac

    output_bam="$output_dir/${filename%.bam}_RGfixed.bam"

    echo "Processing $filename -> $output_bam with sample name $sample"

    # Add read groups using samtools
    samtools addreplacerg -r "@RG\tID:$sample\tSM:$sample\tPL:illumina\tLB:lib1" -o "$output_bam" "$bamfile"

    # Index the new BAM file
    samtools index "$output_bam"
done

echo "✅ Read group fixing complete!"
