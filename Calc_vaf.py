###
# Calculating Variant Allele Frequency from VCF Format
# Author: Rob Schultz
# Purpose: Interview with Leidos Biomedical
###

import os

#check if the output file exists and create it if it does not
output_name = 'Schultz-VAF_Val_ClinOmics_Compass_026_T1D_E.strelka.snvs.raw'
extension = '.vcf'
new_file = output_name + extension
count = 0
while new_file in os.listdir():
    new_file = output_name + str(count) + extension
    count += 1
else:
    with open(new_file,'w') as file:
        pass # Create an empty new file to append to later


isInfoTagAdded = False
vaf_info = '##INFO=<ID=VAF_1,Number=A,Type=Integer,Description=\"Variant Allele Frequency for first sample\">\n##INFO=<ID=VAF_2,Number=A,Type=Integer,Description=\"Variant Allele Frequency for second sample\">\n'
output1 = [] # list of each line from updated VCF for debugging purposes
# Read in the VCF text
with open("Val_ClinOmics_Compass_026_T1D_E.strelka.snvs.raw.vcf",'r') as vcf:
    for line in vcf:
        #line = vcf.readline()
        if line[0] == '#': # if the line is meta-information or header line, do not change it.
            new_line = line
        else: # if the line is a variant line, calculate VAF
            #parse line
            vals = line.split('\t')
            format_col = vals[8].split(':')
            sample_n = vals[9].split(':') # n will notate the normal sample
            sample_t = vals[10].split(':') # t will notate the tumor sample
            format_dict_n = {}
            format_dict_t = {}
            count = 0
            for field in format_col:
                format_dict_n[field] = sample_n[count]
                format_dict_t[field] = sample_t[count]
                count += 1                              #
            #
            ref = vals[3] + 'U'
            v = vals[4].split(',')
            alt = [i for i in map(lambda x: x+'U',v)]
            alt_val_n = [int(format_dict_n[a].split(',')[0]) for a in alt]
            alt_val_t = [int(format_dict_t[a].split(',')[0]) for a in alt]
            # Calculate variant allele frequencies
            ref_val_n = int(format_dict_n[ref].split(',')[0])
            ref_val_t = int(format_dict_t[ref].split(',')[0])
            vaf_n = [i for i in map(lambda s: 0 if s==0 else s/(ref_val_n+s),alt_val_n)]
            vaf_t = [i for i in map(lambda s: 0 if s==0 else s/(ref_val_t+s),alt_val_t)]
            #alt_val_nt = [alt_val_n[i] + alt_val_t[i] for i in range(len(v))]
            #ref_val_nt = ref_val_n +ref_val_t
            #vaf_nt = [i for i in map(lambda s: 0 if s == 0 else s/(ref_val_nt+s),alt_val_nt)] 
            # Update line with the vafs rounded to 2 decimals
            vals[7] = 'VAF_2=' + ','.join([str(round(vaf,2)) for vaf in vaf_t]) + ';' + vals[7]
            vals[7] = 'VAF_1=' + ','.join([str(round(vaf,2)) for vaf in vaf_n]) + ';' + vals[7]
            new_line = '\t'.join(vals)
        with open(new_file,'a') as output:
            # Add the new INFO tags before the first INFO tag
            if line[:6] == "##INFO" and not isInfoTagAdded: 
                output1.append(vaf_info)
                #write the new INFO tag to file
                output.write(vaf_info)
                isInfoTagAdded = True
            output1.append(new_line)
            # Write the updated line to file
            output.write(new_line)
            







    





# Directions from Strelka UserGuide
#Somatic SNVs:
#refCounts = Value of FORMAT column $REF + “U” (e.g. if REF="A" then use the value in FOMRAT/AU)
#altCounts = Value of FORMAT column $ALT + “U” (e.g. if ALT="T" then use the value in FOMRAT/TU)
#tier1RefCounts = First comma-delimited value from $refCounts
#tier1AltCounts = First comma-delimited value from $altCounts
#Somatic allele freqeuncy is $tier1AltCounts / ($tier1AltCounts + $tier1RefCounts)


























