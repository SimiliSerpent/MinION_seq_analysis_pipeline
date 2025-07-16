#!/bin/bash

cat b01_porechopped_map2_multiv3.human.samstats | grep ^SN | cut -f 2- > barcode_01.human.summary
cat b02_porechopped_map2_multiv3.human.samstats | grep ^SN | cut -f 2- > barcode_02.human.summary
cat b03_porechopped_map2_multiv3.human.samstats | grep ^SN | cut -f 2- > barcode_03.human.summary
cat b04_porechopped_map2_multiv3.human.samstats | grep ^SN | cut -f 2- > barcode_04.human.summary
cat b05_porechopped_map2_multiv3.human.samstats | grep ^SN | cut -f 2- > barcode_05.human.summary
cat b06_porechopped_map2_multiv3.human.samstats | grep ^SN | cut -f 2- > barcode_06.human.summary
cat b07_porechopped_map2_multiv3.human.samstats | grep ^SN | cut -f 2- > barcode_07.human.summary
cat unc_porechopped_map2_multiv3.human.samstats | grep ^SN | cut -f 2- > unclassified.human.summary

cat b01_porechopped_map2_multiv3.mouse.samstats | grep ^SN | cut -f 2- > barcode_01.mouse.summary
cat b02_porechopped_map2_multiv3.mouse.samstats | grep ^SN | cut -f 2- > barcode_02.mouse.summary
cat b03_porechopped_map2_multiv3.mouse.samstats | grep ^SN | cut -f 2- > barcode_03.mouse.summary
cat b04_porechopped_map2_multiv3.mouse.samstats | grep ^SN | cut -f 2- > barcode_04.mouse.summary
cat b05_porechopped_map2_multiv3.mouse.samstats | grep ^SN | cut -f 2- > barcode_05.mouse.summary
cat b06_porechopped_map2_multiv3.mouse.samstats | grep ^SN | cut -f 2- > barcode_06.mouse.summary
cat b07_porechopped_map2_multiv3.mouse.samstats | grep ^SN | cut -f 2- > barcode_07.mouse.summary
cat unc_porechopped_map2_multiv3.mouse.samstats | grep ^SN | cut -f 2- > unclassified.mouse.summary

cat b01_porechopped_map2_multiv3.lambda.samstats | grep ^SN | cut -f 2- > barcode_01.lambda.summary
cat b02_porechopped_map2_multiv3.lambda.samstats | grep ^SN | cut -f 2- > barcode_02.lambda.summary
cat b03_porechopped_map2_multiv3.lambda.samstats | grep ^SN | cut -f 2- > barcode_03.lambda.summary
cat b04_porechopped_map2_multiv3.lambda.samstats | grep ^SN | cut -f 2- > barcode_04.lambda.summary
cat b05_porechopped_map2_multiv3.lambda.samstats | grep ^SN | cut -f 2- > barcode_05.lambda.summary
cat b06_porechopped_map2_multiv3.lambda.samstats | grep ^SN | cut -f 2- > barcode_06.lambda.summary
cat b07_porechopped_map2_multiv3.lambda.samstats | grep ^SN | cut -f 2- > barcode_07.lambda.summary
cat unc_porechopped_map2_multiv3.lambda.samstats | grep ^SN | cut -f 2- > unclassified.lambda.summary

cat b01_porechopped_map2_multiv3.SARSCoV2.samstats | grep ^SN | cut -f 2- > barcode_01.SARSCoV2.summary
cat b02_porechopped_map2_multiv3.SARSCoV2.samstats | grep ^SN | cut -f 2- > barcode_02.SARSCoV2.summary
cat b03_porechopped_map2_multiv3.SARSCoV2.samstats | grep ^SN | cut -f 2- > barcode_03.SARSCoV2.summary
cat b04_porechopped_map2_multiv3.SARSCoV2.samstats | grep ^SN | cut -f 2- > barcode_04.SARSCoV2.summary
cat b05_porechopped_map2_multiv3.SARSCoV2.samstats | grep ^SN | cut -f 2- > barcode_05.SARSCoV2.summary
cat b06_porechopped_map2_multiv3.SARSCoV2.samstats | grep ^SN | cut -f 2- > barcode_06.SARSCoV2.summary
cat b07_porechopped_map2_multiv3.SARSCoV2.samstats | grep ^SN | cut -f 2- > barcode_07.SARSCoV2.summary
cat unc_porechopped_map2_multiv3.SARSCoV2.samstats | grep ^SN | cut -f 2- > unclassified.SARSCoV2.summary

cat b01_porechopped_map2_multiv3.Ecoli.samstats | grep ^SN | cut -f 2- > barcode_01.Ecoli.summary
cat b02_porechopped_map2_multiv3.Ecoli.samstats | grep ^SN | cut -f 2- > barcode_02.Ecoli.summary
cat b03_porechopped_map2_multiv3.Ecoli.samstats | grep ^SN | cut -f 2- > barcode_03.Ecoli.summary
cat b04_porechopped_map2_multiv3.Ecoli.samstats | grep ^SN | cut -f 2- > barcode_04.Ecoli.summary
cat b05_porechopped_map2_multiv3.Ecoli.samstats | grep ^SN | cut -f 2- > barcode_05.Ecoli.summary
cat b06_porechopped_map2_multiv3.Ecoli.samstats | grep ^SN | cut -f 2- > barcode_06.Ecoli.summary
cat b07_porechopped_map2_multiv3.Ecoli.samstats | grep ^SN | cut -f 2- > barcode_07.Ecoli.summary
cat unc_porechopped_map2_multiv3.Ecoli.samstats | grep ^SN | cut -f 2- > unclassified.Ecoli.summary

cat b01_porechopped_map2_multiv3.Saureus.samstats | grep ^SN | cut -f 2- > barcode_01.Saureus.summary
cat b02_porechopped_map2_multiv3.Saureus.samstats | grep ^SN | cut -f 2- > barcode_02.Saureus.summary
cat b03_porechopped_map2_multiv3.Saureus.samstats | grep ^SN | cut -f 2- > barcode_03.Saureus.summary
cat b04_porechopped_map2_multiv3.Saureus.samstats | grep ^SN | cut -f 2- > barcode_04.Saureus.summary
cat b05_porechopped_map2_multiv3.Saureus.samstats | grep ^SN | cut -f 2- > barcode_05.Saureus.summary
cat b06_porechopped_map2_multiv3.Saureus.samstats | grep ^SN | cut -f 2- > barcode_06.Saureus.summary
cat b07_porechopped_map2_multiv3.Saureus.samstats | grep ^SN | cut -f 2- > barcode_07.Saureus.summary
cat unc_porechopped_map2_multiv3.Saureus.samstats | grep ^SN | cut -f 2- > unclassified.Saureus.summary
