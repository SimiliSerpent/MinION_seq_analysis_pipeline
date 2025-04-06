# Import required python modules
import os, shutil
import json
import glob
import sys
import argparse
from pathlib import Path
from subprocess import Popen

# Import auxiliary functions
from utils import *

# Abs path to taxonomy file
tax_file='/env/cnrgh/proj/math_stats/scratch/vacus/Projet_CAPASVIR/scratch/Data/taxonomy.tab'

def getTaxonomy(file):
  taxo={}
  with open(file) as f :
    for line in f.readlines():
      taxo[line.split('\t')[0]]= [line.split('\t')[2],line.split('\t')[4].replace('\n','')]
      
  return taxo

def convert(kraken_report, report_krona):
  print('importing taxonomy')
  taxonomy=getTaxonomy(tax_file)
  print('converting data')
  with open(kraken_report, 'r') as f_in :
    lines = f_in.readlines()

  with open(report_krona, 'w') as f_out :
    for line in lines :
        
      if(int(line.split('\t')[2])!=0):
        tax_id=line.split('\t')[4]
        read_max = line.split('\t')[1]
        read_num = line.split('\t')[2]
        init_val = line.split('\t')[-1]
        temp_val = init_val.split('  ')[-1]
        if(tax_id=='0'):
          f_out.write(read_num + '\t' + temp_val)
        else:
          try:
            lignage = [temp_val]
            tmp_tax_id=tax_id
            while tmp_tax_id != '1' :
              parent = taxonomy[tmp_tax_id]
              lignage.append(parent[1])
              tmp_tax_id = parent[0]
            f_out.write(read_num+'\troot')
            if(tax_id=='1'):
              f_out.write('\n')
            else:
              for i in range(len(lignage) -1):
                f_out.write('\t'+lignage[len(lignage)-1 - i])
              f_out.write('\n')
          except:
            print("problem with tax id : "+str(tmp_tax_id)+" "+temp_val+" probably not in taxonomy")
  # os.chmod(report_krona,0o777)
  print("done converting")
  return 0


def convert_to_krona(report, outdir):
  print('converting:', report, 'to', outdir)

  if not os.path.exists(report + ".html"):

    krona_file=report.replace(":", "_") + ".krona"
    convert(report, krona_file)

    cmd = [
      '/env/cnrgh/proj/math_stats/scratch/vacus/Projet_CAPASVIR/scratch/Analysis/scripts/KronaTools-2.8.1/bin/ktImportText',
      os.path.abspath(krona_file),
      '-o',
      os.path.abspath(outdir.rstrip('/') + '/' + report.split('/')[-1] + '.html')
    ]
    print(cmd)
    p = Popen(cmd)
    p.communicate()
    
    os.remove(krona_file)


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-r', '--reports_dir', type=str, required=True,
                      help='Path to the kraken reports directory')
  parser.add_argument('-o', '--output_dir', type=str, default='.',
                      help='Path to the output file directory'
                      + ' (default: current/working directory)')

  args = parser.parse_args()

  # Retrieve kraken report files to convert
  report_files = find_file(args.reports_dir, '.k2_report')
  report_files = [args.reports_dir.rstrip('/') + ('/') + file for file in report_files]
  print('report files found:', report_files)

  for file in report_files:
    convert_to_krona(file, args.output_dir)

  return 0


if __name__ == "__main__":
  sys.exit(main())
