# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 16:49:46 2021

@author: dmitr
"""

import re
import pandas as pd
import openpyxl
import os
import sys
import time

class Molcas:
    def __init__(self, molfile, path):
        if path[-1] != os.sep:
            self.PATH = path + os.sep
        else:
            self.PATH = path
        if re.findall(r'.log', molfile):
            self.MOL = self.PATH + molfile
        else:
            self.MOL = self.PATH + molfile + '.log'

        self.Excel = [''.join((self.MOL.replace('.log',''), '_RASSCF.xlsx'))]
        print(self.Excel)

        self.RASnum = 0
        self.CIroots = []
        # Number of inactive and active orbitals
        self.numIAO = []
        self.coords = {}

    def сreate_Excel(self):
        for Excel_file in self.Excel:
            xlsxfile = pd.ExcelWriter(Excel_file, engine='xlsxwriter')
            xlsxfile.close()
        return print('Excel file was created')
    
    def write_to_Excel(self, Excel_file, sheet_name, *args):
        for data_frame in args:
            with pd.ExcelWriter(Excel_file, engine='openpyxl', mode='a') as writer:
                data_frame.to_excel(writer, sheet_name= sheet_name)
        return 'DataFrame is written to Excel file'
    
    def delete_sheet_Excel(self):
        for Excel_file in self.Excel:
            xlsxfile = openpyxl.load_workbook(Excel_file)
            xlsxfile.remove(xlsxfile['Sheet1'])
            xlsxfile.save(Excel_file)
            xlsxfile.close()
        return 'Sheet was deleted'

    #Check of the convergence of RASSCF block
    def convergence_check(self):
        conv = False
        with open(self.MOL, 'r') as MOL:
            for line in MOL:
                if not re.search(r'Convergence after', line):
                    continue
                else:
                    conv = True
                    print('Calculation converged.')
                    return conv
            if not conv:
                print('No convergence.')
                return conv

    # Number of RASSCF blocks and CI roots for each block
    def get_RASSCFandCI_numbers(self):
        with open(self.MOL,'r') as MOL:
            CI = False
            for line in MOL:
                if re.search(r'&RASSCF', line):
                    self.RASnum += 1
                    continue
                if re.search(r'CIROOT',line) and len(line.split()) == 1:
                    CI = True
                    continue
                if re.search(r'CIROOT',line) and len(line.split()) > 1:
                    self.CIroots.append(int(line.split()[1]))
                    continue
                if CI:
                    self.CIroots.append(int(line.split()[0]))
                    CI = False
                    continue
                if re.search(r'(\(\)){6}', line):
                    break
        return [self.RASnum, self.CIroots]
    
    def extract_geometry(self):
        with open(self.MOL, 'r') as MOL:
            coords = {'Atoms': [], 'X': [], 'Y': [], 'Z': []}
            g_block = False
            for line in MOL:
                if re.search(r'Label {8}X {12}Y {12}Z ', line):
                    g_block = True
                    continue
                if g_block and line.split()[0].isnumeric():
                    coords['Atoms'].append(line.split()[1])
                    coords['X'].append(float(line.split()[2]))
                    coords['Y'].append(float(line.split()[3]))
                    coords['Z'].append(float(line.split()[4]))
                    continue
                if g_block and re.search('Nuclear repulsion energy', line):
                    break
            self.coords = coords
            del coords
        return print('Geometry is stored')

    # Division of the molcas output into pieces with different RASSCF blocks
    def divide_RASSCF_output(self):
        i = 1
        with open(self.MOL.replace('.log','_'+str(i)+'.log'), 'w') as MOL1:
            MOL1.write('')
        with open(self.MOL, 'r') as MOL:
            for line in MOL:
                if i == 1 and MOL1.closed:
                    MOL1 = open(self.MOL.replace('.log','_'+str(i)+'.log'), 'w')
                if i%2 != 0:
                    if not (re.search(r'CMON in NATORB_RASSCF after ORDER_ARRAYS', line)
                            or re.search(r'Final results', line)):
                        MOL1.write(line)
                        continue
                    if re.search(r'Final results', line):
                        MOL1.close()
                        i += 1
                        MOL1 = open(self.MOL.replace('.log','_'+str(i)+'.log'), 'w')
                        MOL1.write(line)
                        continue
                    if re.search(r'CMON in NATORB_RASSCF after ORDER_ARRAYS', line):
                        MOL1.close()
                        i += 1
                        continue
                else:
                    if re.search(r'Final results', line) and MOL1.closed:
                        MOL1 = open(self.MOL.replace('.log','_'+str(i)+'.log'), 'w')
                        MOL1.write(line)
                        continue
                    if not re.search(r'Module rasscf spent',line) and not MOL1.closed:
                        MOL1.write(line)
                        continue
                    if re.search(r'Module rasscf spent',line) and not MOL1.closed:
                        MOL1.write(line)
                        MOL1.close()
                        if i == 2*self.RASnum:
                            break
                        else:
                            i += 1
                            MOL1 = open(self.MOL.replace('.log','_'+str(i)+'.log'), 'w')
                            continue
        return print(self.MOL + ' is divided')
    
    # Writing to tables CI coefficients of every root in RASSCF blocks
    def get_CIroots(self):
        roots = pd.DataFrame({})
        CIs = {'Energy': [], 'Orbitals': [], 'Coeff': [], 'Weight': []}
        for i in range(self.RASnum):
            RasEnergies = {'RasEnergies': []}
            MOL = open(self.MOL.replace('.log','_'+str(2*i+2)+'.log'),'r')
            for line in MOL:
                if re.search(r':: {4}RASSCF root number',line):
                    RasEnergies['RasEnergies'].append(line.split()[-1])
                    continue
                if re.search(r'\++ {4}Molecular orbitals',line):
                    roots = pd.concat([roots, pd.DataFrame(RasEnergies)], axis=1, join='outer')
                    break
            MOL.close()
            MOL = open(self.MOL.replace('.log','_'+str(2*i+1)+'.log'),'r')
            ROOT = False
            for line in MOL:
                if re.search(r'energy=',line):
                    CIs['Energy'] = [float(line.split()[1])]
                    ROOT = True
                    continue
                if ROOT and line.split() and line.split()[0].isnumeric():
                    CIs['Orbitals'].append(line.split()[1])
                    CIs['Coeff'].append(float(line.split()[2]))
                    CIs['Weight'].append(float(line.split()[3]))
                    continue
                if ROOT and not line.split():
                    ROOT = False
                    CIs['Energy'] += [None for i in range(1,len(CIs['Coeff']))]
                    roots = pd.concat([roots, pd.DataFrame(CIs)], axis=1, join='outer')
                    CIs = {'Energy': [], 'Orbitals': [], 'Coeff': [], 'Weight': []}
                    continue
                if re.search(r'Natural orbitals and occupation numbers for root', line):
                    self.write_to_Excel(self.Excel[0], 'RASSCF_module_'+str(i+1), roots)
                    roots = pd.DataFrame({})
                    break
            MOL.close()
        return print('CI data is written')

    def get_orbs_number(self):
        with open(self.MOL.replace('.log','_2.log')) as MOL1: #Get number of inactive and active orbitals
            for line in MOL1:
                if re.search(r'Number of inactive orbitals', line):
                    self.numIAO.append(int(line.split()[4]))
                    continue
                elif re.search(r'Number of active orbitals', line):
                    self.numIAO.append(int(line.split()[4]))
                    break
                else: continue
        return print('Number of inactive and active orbitals is retrieved')

    def get_active_orbs(self, active_orbs: list, num: int):
        with open(self.MOL.replace('.log', '_' + str(2 * (num + 1)) + '.log'), 'r') as MOL1:
            j = 0
            for line in MOL1:
                if j == self.numIAO[1] and len(line.split()) == 3:
                    del j
                    break
                if line.split() and j != 0 and line.split()[1] != '0.0000':
                    active_orbs[j - 1].append(line.replace('(', ' ').replace(')', ' '))
                    continue
                if line.split() and line.split()[0] == str(self.numIAO[0] + j + 1) and line.split()[1] == '0.0000':
                    active_orbs[j].append(float(line.split()[2]))  # add a number of electrons at orbital to list
                    j += 1
                    continue
        return active_orbs

    def get_active_space(self):
        self.get_orbs_number()

        for i in range(self.RASnum):
            active_orbs = [[str(self.numIAO[0] + j + 1)] for j in range(self.numIAO[1])]
            active_space_dict = {}
            active_space_pd = pd.DataFrame(active_space_dict)

            self.get_active_orbs(active_orbs, i)

            for j in range(len(active_orbs)): #filling the dictionary with active space data with following write to dataframe
                active_space_dict[active_orbs[0][0]] = [active_orbs[0][1]] #dictionary first kv: orbital-number of electrons
                active_space_dict['Atom'] = []
                active_space_dict['Orbital'] = []
                active_space_dict['WFcoeff'] = []
                active_space_dict['WFCin2'] = []
                for k in range(2,len(active_orbs[0])):
                    line = active_orbs[0][k].split()
                    active_space_dict['Atom'] += line[1:len(line):4]
                    active_space_dict['Orbital'] += line[2:len(line):4]
                    active_space_dict['WFcoeff'] += [float(l) for l in line[3:len(line):4]]
                active_space_dict['WFCin2'] += [l**2 for l in active_space_dict['WFcoeff']]
                active_space_dict['WFCin2'] = [l/sum(active_space_dict['WFCin2'],0) for l in active_space_dict['WFCin2']]
                active_space_dict[active_orbs[0][0]] += [None for l in range(1,len(active_space_dict['Atom']))]
                active_space_pd = pd.concat([active_space_pd, pd.DataFrame(active_space_dict)], axis=1, join = 'outer')
                del(active_orbs[0])
                active_space_dict = {}
            self.write_to_Excel(self.Excel[0], 'RASSCF_AS_'+str(i+1), active_space_pd)
        return print('Atomic orbitals contributions to AS are written')

    def get_Mulliken_LoProp(self): #Get Mulliken and LoProp properties
        for i in range(self.RASnum):

            mulo_dict = ['Mulliken', 'Mulliken spin', 'LoProp']
            dicts = ['MulCharges', 'MulSpin', 'LPCharges']
            search_patterns = [['N-E', 'Mulliken Bond Order analysis'],
                               ['Total', 'Total electronic spin='],
                               ['Total', 'Natural Bond Order analysis for root']]

            with open(self.MOL.replace('.log','_'+str(2*(i+1))+'.log'), 'r') as MOL1:

                cond_check = False
                num = 0
                root = 1
                population_dict = {'Atoms': self.coords['Atoms']}

                for line in MOL1:
                    if root == self.CIroots[i] + 1:
                        break
                    if re.search(fr"{mulo_dict[num]}" + r" population analysis for root", line, re.IGNORECASE):
                        cond_check = True
                        population_dict[f"{dicts[num]}" + " of root " + str(root)] = []
                        continue
                    if cond_check and re.search(fr"{search_patterns[num][0]}", line, re.IGNORECASE) \
                            and not re.search(r'Total electronic spin=', line, re.IGNORECASE):
                        population_dict[f"{dicts[num]}" + " of root " + str(root)] += line.split()[1:]
                        continue
                    if cond_check and re.search(fr"{search_patterns[num][1]}", line, re.IGNORECASE):
                        population_dict[f"{dicts[num]}" + " of root " + str(root)] = \
                            list(map(float, population_dict[f"{dicts[num]}" + " of root " + str(root)]))
                        cond_check = False
                        if num == 2:
                            num = 0
                            root += 1
                        else:
                            num += 1
                        continue

                self.write_to_Excel(self.Excel[0],'Population_block_'+str(i+1), pd.DataFrame(population_dict))

        return print('Mulliken and LoProp properties are written')

def main():
    start_time = time.time()

    try:
        molcas_file_name, path_file = sys.argv[1], os.getcwd()
    except:
        print("Error! There is no file name.")
    else:
        mf = Molcas(molcas_file_name, path_file)

        if mf.convergence_check():
            mf.сreate_Excel()
            mf.extract_geometry()
            print(mf.get_RASSCFandCI_numbers())
            mf.divide_RASSCF_output()
            mf.get_CIroots()
            mf.get_active_space()
            mf.get_Mulliken_LoProp()
            mf.delete_sheet_Excel()

    time_proc = time.time() - start_time
    print('time spent =', time_proc)

if __name__ == '__main__':
    main()
