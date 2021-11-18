#!/usr/bin/python
#
# print a summary HTML page for all files in result directory
#
# Usage: PROG DIR
# e.g. ~/projects/Larpack/scripts/summaryHTML.py . > summary.html
#
# Mario Fasold, 2011/06/29
from __future__ import division # no surprises for division
import sys
import os
import csv

# Collect list of valid hook directories
base_dir = sys.argv[1]
hook_dirs = [name for name in os.listdir(base_dir) 
             if os.path.isdir(os.path.join(base_dir, name)) and os.path.exists(os.path.join(base_dir, name, "dataLogger.log"))]

hook_dirs = sorted(hook_dirs)


def read_param_value_pair_file(filename, delimiter='\t'):
    pair_dict = {}
    info_file = csv.reader(open(filename, "rU"), delimiter=delimiter)
    for row in info_file:
        if len(row) > 1:
            pair_dict[row[0]] = row[1].strip()
    return pair_dict


# Print HTML header 
print '''<html> 
<head> <title>Hook Summary</title> </head> 
<body bgcolor=#FFFFFF >
<H1 ALIGN=CENTER >Hook Summary</H1>
<style type=text/css> 
p{ margin-top: 1px; margin-bottom: 1px; padding-left: 10px; text-indent: -10px } 
table { text-align:center; table-layout:fixed}
</style> '''


# this is for the most simple pagination there is....
if len(sys.argv) > 2:
    pages_to_print = [int(sys.argv[2])]
else:
    pages_to_print = [1,2,3]



if 1 in pages_to_print:

    for index, hook_dir in enumerate(hook_dirs):
        print "<h2>Chip", index + 1, ": ", hook_dir, "</h2>"
        print "<center><table border=1><tr><td>HookCurve</td><td>HookCurve NS Corrected</td></tr>"
        print '<tr><td><a href="' + hook_dir + '/graphics/Hookplot-Primary.png"><img src="' + hook_dir + '/graphics/Hookplot-Primary.png" border=0 width=400 height=300 ></a></td>'
        print '<td><a href="' + hook_dir + '/graphics/Hookplot-Corrected.png"><img src="' + hook_dir + '/graphics/Hookplot-Corrected.png" border=0 width=400 height=300 ></a></td>'
        print '</tr></table></center>'    


    print '''<h2>Table 1</h2>
<center><table border=1><tr><td></td><td>lg(oBG)</td><td>%N</td><td>lg(N)</td><td>lg(M)</td><td>Alpha</td><td>av R</td></tr>'''

    for index, hook_dir in enumerate(hook_dirs):
        hook_params = read_param_value_pair_file(os.path.join(base_dir, hook_dir, "dataLogger.log"))
        print '<tr><td>Chip ', index + 1,': ', hook_dir, '</td>'
        print '<td>', hook_params["backgroundSubtractionBgMean"], '</td>'
        print '<td>', hook_params["probesetNsFractionCorrected"], '</td>'
        print '<td>', hook_params["hookcurveKinkPointCorrected"], '</td>'
        print '<td>', hook_params["saturationImax"], '</td>'
        print '<td>', hook_params["saturationA"], '</td>'
        print '<td>', hook_params["rGreaterThan05Mean"], '</td>'
    #    print '<td>', hook_params[""], '</td>'
        print '</tr>'

    print '''</table></center> 
<br><center><table border=1>
<tr><td>lg(oBG)</td><td>optical background in log-scale (Affy zone algorithm)</td></tr>
<tr><td>%N</td><td>percent "absent" probes which are predominantly hybridized non-specifically</td></tr>
<tr><td>lg(N)</td><td>mean intensity of the N-background (log scale)</td></tr>
<tr><td>lg(M)</td><td>mean saturation intensity (log scale)</td></tr>

<tr><td>Alpha</td><td>PM/MM-gain in the S-range (asymptotic "height" of the hook)</td></tr>
<tr><td>av R</td><td>mean S/N-ratio</td></tr>
</table></center><br>
'''



    print '''<h2>Table 2</h2>
<center><table border=1><tr><td></td><td>lg(n)</td><td>SD(N)</td><td>av lg(R)</td><td>Beta</td><td>Fi</td></tr>'''

    for index, hook_dir in enumerate(hook_dirs):
        hook_params = read_param_value_pair_file(os.path.join(base_dir, hook_dir, "dataLogger.log"))
        print '<tr><td>Chip ', index + 1,': ', hook_dir, '</td>'
        print '<td>', hook_params["nsMeanProbessetsDiffCorrected(logB)"], '</td>'
        print '<td>', hook_params["hookcurveNsRangeWidthUncorrected"], '</td>'
        print '<td>', hook_params["rGreaterThan05LogPlus1Mean"], '</td>'
        print '<td>', hook_params["saturationF"], '</td>'
        print '<td>', "?", '</td>'
    #    print '<td>', hook_params[""], '</td>'
        print '</tr>'

    print '''</table></center> 
<br><center><table border=1>
<tr><td>lg(n)</td><td>PM/MM-gain in the N-range</td></tr>
<tr><td>SD(N)</td><td>width of the N-range</td></tr>

<tr><td>av lg(R)</td><td>logarithmic mean of the S/N-ratio (Av lg(R+1) )</td></tr>
<tr><td>Beta</td><td>N-binding strength ("width" of the hook) (log scale)</td></tr>
<tr><td>Fi</td><td>mean expression (log scale)</td></tr>
</table></center><br>'''
# print '''
# <h2>Table 3</h2>
# <center><table border=1><tr><td></td><TD colspan="3">PM uncorrected</td><TD colspan="3">PM corrected</td><td>D</td></tr>
# <tr><td></td><td>&mu; Probesets</td><td>&sigma; Probes</td><td>&rho; NS Probes</td><td>&mu; Probesets</td><td>&sigma; Probes</td><td>&rho; NS Probes</td><td></td></tr>'''



if 2 in pages_to_print:
    for index, hook_dir in enumerate(hook_dirs):
        print '<h2>Chip', index + 1, ': ', hook_dir, '</h2>'
        print '<center><table border=1><tr><td>PM-NS</td><td>MM-NS</td><td>PM-S</td><td>MM-S</td></tr>'
        print '<tr><td><img src="' + hook_dir +'/graphics/SensitivityProfileNsPm-Corrected.png" border=0 width=400 height=300 ></td>'
        print '<td><img src="' + hook_dir +'/graphics/SensitivityProfileNsMm-Corrected.png" border=0 width=400 height=300 ></td>'
        print '<td><img src="' + hook_dir +'/graphics/SensitivityProfileSPm-Corrected.png" border=0 width=400 height=300 ></td>'
        print '<td><img src="' + hook_dir +'/graphics/SensitivityProfileSMm-Corrected.png" border=0 width=400 height=300 ></td>'
        print '</tr></table></center>'

    # print '<center><table border=1><tr><td>PM-NS</td><td>MM-NS</td><td>PM-S</td><td>MM-S</td></tr>'    
    # for index, hook_dir in enumerate(hook_dirs):
    #     #print '<h2>Chip', index + 1, ': ', hook_dir, '</h2>'
    #     print '<tr><td><img src="' + hook_dir +'/graphics/SensitivityProfileNsPm-Corrected.png" border=0 width=400 height=300 ></td>'
    #     print '<td><img src="' + hook_dir +'/graphics/SensitivityProfileNsMm-Corrected.png" border=0 width=400 height=300 ></td>'
    #     print '<td><img src="' + hook_dir +'/graphics/SensitivityProfileSPm-Corrected.png" border=0 width=400 height=300 ></td>'
    #     print '<td><img src="' + hook_dir +'/graphics/SensitivityProfileSMm-Corrected.png" border=0 width=400 height=300 ></td>'
    #     print '</tr'
    # print '</table></center>'


    # TODO: add gstacks variant

# Create 3' 5' Images

if 3 in pages_to_print:
    for index, hook_dir in enumerate(hook_dirs):
        print '<h2>Chip', index + 1, ': ', hook_dir, '</h2>'
        print '<center><table border=1><tr><td>Polymerase Bias</td></tr>'

        print '<tr><td><img src="' + hook_dir +'/graphics/PolymeraseBias.png" border=0 width=640 height=480 ></td></tr></table></center>'


print '</body></html>'


