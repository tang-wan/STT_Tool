import os
import numpy as np # type: ignore
import re
from colorama import Fore, Back # type: ignore
import matplotlib.pyplot as plt
# import matplotlib.animation as ama
# from tqdm.autonotebook import tqdm
from itertools import combinations_with_replacement

# =========================
def MKdir(path, rm):
    if os.path.isdir(path) and rm:
        os.system(f"rm -r {path}")
        os.mkdir(path)
    elif os.path.isdir(path) and not(rm):
        pass
    else:
        os.mkdir(path)
# =========================
def Check_out_Word(word: str):
    print(Fore.BLACK, Back.RED + word + Back.RESET, Fore.RESET)
def Process_Word(word: str):
    print(Fore.YELLOW, word , Fore.RESET)
# =========================
def M_Rotate(deg):
    rad = deg/180*np.pi
    return np.array([[np.cos(rad), -np.sin(rad), 0 ],
                     [np.sin(rad),  np.cos(rad), 0 ],
                     [     0     ,       0     , 1 ]])

class POSCARConvert():
    def __init__(self, File, atom_type):
        self.File = File
        self.atom_type = atom_type

        # atom type in str
        self.atom_st          = re.sub(u'([^\u0041-\u007a])', ' ', atom_type)
        # atom type in lsit
        self.atom             = self.atom_st.split()
        
        # atom number in str
        self.atom_num_st      = re.sub(u'([^\u0030-\u0039])', ' ', atom_type).replace(" ", "", 1)
        # atom number in list
        self.atom_num         = list(map(int, self.atom_num_st.split()))
        self.m                = sum(self.atom_num)
        
    def ReadPOSCAR(self):
        File_Name   = self.File
        atom_st     = self.atom_st 
        atom        = self.atom    
        atom_num_st = self.atom_num_st 
        atom_num    = self.atom_num    
        m           = self.m
    # --------------
        with open(File_Name, 'r') as R_File:
            R_content = R_File.readlines()
        L_const = float(R_content[1].split()[0])
        a_vec   = np.float64(R_content[2].split())*L_const
        b_vec   = np.float64(R_content[3].split())*L_const
        c_vec   = np.float64(R_content[4].split())*L_const
        atom_pos_list = R_content[8:8+m]
        LatticeConstant = (np.sum(a_vec[0:2]**2))**0.5

    # --------------
        ## atom position in direction
        atom_array_dir = np.zeros(3)
        for pos in atom_pos_list:
            atom_dir = np.float64(pos.split())
            atom_array_dir = np.vstack((atom_array_dir,
                                        atom_dir)
                                        )
        atom_array_dir = np.delete(atom_array_dir, 0, axis=0)
    # --------------
        ## atom position in Cartition
        atom_array_car = np.zeros(3)
        for a in atom_array_dir:

            atom_array_car = np.vstack((atom_array_car, 
                                        np.round(a[0]*a_vec, 15) + np.round(a[1]*b_vec, 15) + np.round(a[2]*c_vec, 15))
                                        )
        atom_array_car = np.delete(atom_array_car, 0, axis=0)
    # --------------
        ## Making atom list
        atom_list = []
        for repeat in range(len(atom_num)):
            repeat_num = atom_num[repeat]
            for i in range(repeat_num):
                atom_list.append(atom[repeat])
        atom_array = np.array(atom_list)
    # ============== #
        self.atom_array = atom_array
        self.a_vec, self.b_vec, self.c_vec = a_vec, b_vec, c_vec
        self.LatticeConstant = LatticeConstant
        self.atom_array_dir  = atom_array_dir
        self.atom_array_car  = atom_array_car

        return atom_array, (a_vec, b_vec, c_vec), LatticeConstant, atom_array_dir, atom_array_car
    
    def PrintPOSCAR(self):
        print(self.atom_array)
        print(self.LatticeConstant)
        print(self.a_vec, self.b_vec, self.c_vec)
        print()
        print('Direct')
        print(self.atom_array_dir)
        print()
        print('Cart')
        print(self.atom_array_car)
# ----------    
    def atoms_xyz(self, r, theta, phi):
        atom_array = self.atom_array
        atom_array_car = self.atom_array_car
        
        r_list = []
        for mo, sub_r in np.c_[self.atom_num, r]:
            # print(mo, sub_r)
            for _ in range(int(mo)):
                r_list.append(float(sub_r))
        r_array = np.array(r_list)
                

        with open('atoms.xyz', 'w') as W_File:
            lines = [
               f'{len(self.atom_array_car)}\n'
               f"AtomType	{'X':^18}	{'Y':^18}	{'Z':^18}	  SpinPolarization_r  SpinPolarization_theta  SpinPolarization_phi\n"
            ]
            W_File.writelines(lines)

        with open('atoms.xyz', 'a') as A_File:
            for atom, pos_x, pos_y, pos_z, r in np.c_[atom_array, atom_array_car, r_array]:
                lines = [
                    f'{atom:^8}   {pos_x:<018}   {pos_y:<018}   {pos_z:<018}   {float(r):<018}    {float(theta):<018}     {float(phi):<018}\n'
                ]
                A_File.writelines(lines)

    def BtoTF(self, repeatNum=1, vacLen=10):
        
        a_vec, b_vec, c_vec = self.a_vec, self.b_vec, self.c_vec
        atom_array_dir = self.atom_array_dir
        height = c_vec[-1]
        Zscale = ((repeatNum)*height+vacLen)/height

        self.new_a_vec = a_vec
        self.new_b_vec = b_vec
        self.new_c_vec = c_vec*(1, 1, Zscale)

        new_atom_array_dir = atom_array_dir
        for i in range(repeatNum):
            shift_array = np.array([[0, 0, i], 
                                    [0, 0, i],
                                    [0, 0, i]])
            shifted_array = atom_array_dir + shift_array
            if i == 0:
                pass
            else:
                new_atom_array_dir = np.append(new_atom_array_dir, shifted_array)
                new_atom_array_dir = new_atom_array_dir.reshape(3*(i+1), 3)
                # print(new_atom_array_dir)

        new_atom_array_dir = new_atom_array_dir*np.array([[1, 1, 1/Zscale]]*(3*repeatNum))
        shift = 0.5-np.mean(new_atom_array_dir[:,2])
        self.new_atom_array_dir = new_atom_array_dir + np.array([[0, 0, shift]]*(3*repeatNum))

# ----------
    def Rotate_deg(self, deg, Print):
        a_vec, b_vec, c_vec = self.a_vec, self.b_vec, self.c_vec
        atom_array_dir = self.atom_array_dir

        self.new_a_vec = a_vec
        self.new_b_vec = b_vec
        self.new_c_vec = c_vec

        self.new_a_vec = np.round(np.transpose(M_Rotate(deg)@np.transpose(self.new_a_vec)), 14)
        self.new_b_vec = np.round(np.transpose(M_Rotate(deg)@np.transpose(self.new_b_vec)), 14)
        self.new_c_vec = np.round(np.transpose(M_Rotate(deg)@np.transpose(self.new_c_vec)), 14)

        new_atom_array_dir = np.zeros(3)
        for pos in atom_array_dir:
            
            new_atom_array_dir = np.vstack((new_atom_array_dir,
                                            np.array([pos[0], pos[1], pos[2]]))
                                            )

        new_atom_array_dir = np.delete(new_atom_array_dir, 0, axis=0)
        self.new_atom_array_dir = new_atom_array_dir
        
        if Print:
            print()
            print(f"Hexagonal lattice constant, rotated with {deg} degree")
            print(new_atom_array_dir)
        else:
            pass
    
    def RepeatN(self, repeatVec, Print):
        self.repeatVec = repeatVec

        repeatX = repeatVec[0]
        repeatY = repeatVec[1]
        repeatZ = repeatVec[2]

        a_vec, b_vec, c_vec = self.a_vec, self.b_vec, self.c_vec
        atom_array_dir = self.atom_array_dir

        self.new_a_vec = a_vec*repeatX
        self.new_b_vec = b_vec*repeatY
        self.new_c_vec = c_vec*repeatZ

        new_atom_array_dir = np.zeros(3)
        
        for pos in atom_array_dir:
            for x in range(repeatX):
                for y in range(repeatY):
                    for z in range(repeatZ):            
                        new_atom_array_dir = np.vstack((new_atom_array_dir,
                                                        np.array([(pos[0]+x)/repeatX, (pos[1]+y)/repeatY, (pos[2]+z)/repeatZ]))
                                                        )
                        
        new_atom_array_dir = np.delete(new_atom_array_dir, 0, axis=0)        

        self.new_atom_array_dir = new_atom_array_dir

        if Print:
            print()
            print("Tetragonal lattice constant u1")
            print(new_atom_array_dir)
        else:
            pass

# ----------
    def HextoTri(self, Print=True):
        a_vec, b_vec, c_vec = self.a_vec, self.b_vec, self.c_vec
        atom_array_dir = self.atom_array_dir

        self.new_a_vec = a_vec
        self.new_b_vec = np.transpose(M_Rotate(30)@np.transpose(b_vec))*np.sqrt(3)
        self.new_c_vec = c_vec

        self.new_a_vec = np.round(np.transpose(M_Rotate(30)@np.transpose(self.new_a_vec)), 14)
        self.new_b_vec = np.round(np.transpose(M_Rotate(30)@np.transpose(self.new_b_vec)), 14)
        self.new_c_vec = np.round(np.transpose(M_Rotate(30)@np.transpose(self.new_c_vec)), 14)

        new_atom_array_dir = np.zeros(3)
        for pos in atom_array_dir:
            
            if abs(pos[0]-0) < 1e-6:
                # print('A')
                new_atom_array_dir = np.vstack((new_atom_array_dir,
                                                np.array([pos[0], pos[1], pos[2]]))
                                                )
                new_atom_array_dir = np.vstack((new_atom_array_dir,
                                                np.array([0.5, 0.5, pos[2]]))
                                                )
            elif abs(pos[0]-1/3) < 1e-6:
                # print('C')
                new_atom_array_dir = np.vstack((new_atom_array_dir, 
                                                np.array([1/2, 1/6, pos[2]]))
                                                )
                new_atom_array_dir = np.vstack((new_atom_array_dir, 
                                                np.array([0, (1/2+1/6), pos[2]]))
                                                )

            elif abs(pos[0]-2/3) < 1e-6:
                # print('B')
                new_atom_array_dir = np.vstack((new_atom_array_dir, 
                                                np.array([0, 1/3, pos[2]]))
                                                )
                new_atom_array_dir = np.vstack((new_atom_array_dir, 
                                                np.array([0.5, (1/2+1/3), pos[2]]))
                                                )
            else:
                print("This is not HCP")
                break
        new_atom_array_dir = np.delete(new_atom_array_dir, 0, axis=0)
        self.new_atom_array_dir = new_atom_array_dir
        
        if Print:
            print()
            print("Tetragonal lattice constant u1")
            print(new_atom_array_dir)
        else:
            pass

    def WritePOSCAR(self, POSCARname , 
                    new_atom_type, 
                    ):
        
        new_atom_st = re.sub(u'([^\u0041-\u007a])', ' ', new_atom_type)
        new_atom_num_st = re.sub(u'([^\u0030-\u0039])', ' ', new_atom_type).replace(" ", "", 1)
        
        with open(f'{POSCARname}', 'w') as ini_POSCAR:

            a1x, a1y, a1z = self.new_a_vec[0], self.new_a_vec[1], self.new_a_vec[2]
            a2x, a2y, a2z = self.new_b_vec[0], self.new_b_vec[1], self.new_b_vec[2]
            a3x, a3y, a3z = self.new_c_vec[0], self.new_c_vec[1], self.new_c_vec[2]

            lines = [
                f'{new_atom_type}\n',
                f"{' ':<3}{1.0:<016}\n",
                f"{' ':<4}{a1x:<019}{' ':<4}{a1y:<019}{' ':<4}{a1z:<019}\n",
                f"{' ':<4}{a2x:<019}{' ':<4}{a2y:<019}{' ':<4}{a2z:<019}\n",
                f"{' ':<4}{a3x:<019}{' ':<4}{a3y:<019}{' ':<4}{a3z:<019}\n",
                f"{' ':<3}{new_atom_st}\n",
                f"{' ':<3}{new_atom_num_st}\n",
                'Direct\n',
                
            ]
            ini_POSCAR.writelines(lines)
        
        with open(f'{POSCARname}', 'a') as A_File:
            for pos in self.new_atom_array_dir:
                lines=[
                    f"{' ':<2}{pos[0]:<019}{' ':<2}{pos[1]:<019}{' ':<2}{pos[2]:<019}\n",
                       ]
                A_File.writelines(lines)

class POSCARRead():
    def __init__(self, path):
        with open(path, "r") as R_File:
            dataALL = R_File.readlines()
        Const = np.float64(dataALL[1])
        Vector_a1 = np.float64(dataALL[2].split())*Const
        Vector_a2 = np.float64(dataALL[3].split())*Const
        Vector_a3 = np.float64(dataALL[4].split())*Const
        LatticeVector = np.array([Vector_a1, 
                                Vector_a2, 
                                Vector_a3])
        LatticeConstant_a1 = np.sqrt(np.sum(Vector_a1*Vector_a1))
        LatticeConstant_a2 = np.sqrt(np.sum(Vector_a2*Vector_a2))
        LatticeConstant_a3 = np.sqrt(np.sum(Vector_a3*Vector_a3))
        print("Lattice Constant:")
        print(LatticeConstant_a1, LatticeConstant_a2, LatticeConstant_a3)
        print("---")
        print("Lattice Vector:")
        print(LatticeVector)
        print("=====")
        AtomNum = np.int32(dataALL[6].split())
        AtomPos = np.array([])
        for a in range(sum(AtomNum)):
            try:
                AtomPos = np.append(AtomPos, np.float64(dataALL[8+a].split()))
            except ValueError:
                AtomPos = np.append(AtomPos, np.float64(dataALL[8+a].split()[:-1]))
        AtomPos = AtomPos.reshape(sum(AtomNum), 3)
        print("Atom Position:")
        print(AtomPos)
        # ---
        self.Vector_a1, self.Vector_a2, self.Vector_a3 = Vector_a1, Vector_a2, Vector_a3
        self.LatticeConstant_a1, self.LatticeConstant_a2, self.LatticeConstant_a3 = LatticeConstant_a1, LatticeConstant_a2, LatticeConstant_a3
        self.AtomPos = AtomPos

    def _ReadPos(self, TargetAtom_1, TargetAtom_2, p=True):
        Vector_a1, Vector_a2, Vector_a3 = self.Vector_a1, self.Vector_a2, self.Vector_a3
        LatticeConstant_a1, LatticeConstant_a2, LatticeConstant_a3 = self.LatticeConstant_a1, self.LatticeConstant_a2, self.LatticeConstant_a3
        AtomPos = self.AtomPos

        AtomPos_1_Frac = AtomPos[TargetAtom_1]
        AtomPos_2_Frac = AtomPos[TargetAtom_2]
        if p:
            print("Atom1 Position (Fractional):")
            print(AtomPos_1_Frac)
            # print("---")
            print("Atom2 Position:")
            print(AtomPos_2_Frac)
            print("---")
        AtomPos_1_Cart = AtomPos_1_Frac[0]*Vector_a1 + AtomPos_1_Frac[1]*Vector_a2 + AtomPos_1_Frac[2]*Vector_a3
        AtomPos_2_Cart = AtomPos_2_Frac[0]*Vector_a1 + AtomPos_2_Frac[1]*Vector_a2 + AtomPos_2_Frac[2]*Vector_a3
        if p:
            print("Atom1 Position (Cartesian):")
            print(AtomPos_1_Cart)
            # print("---")
            print("Atom2 Position:")
            print(AtomPos_2_Cart)
            print("---")

        return AtomPos_1_Cart, AtomPos_2_Cart

    def LayerHeight(self, TargetAtom_1, TargetAtom_2, p=True): # Atom2 - Atom1

        AtomPos_1_Cart, AtomPos_2_Cart = self._ReadPos(TargetAtom_1, TargetAtom_2, p=p)
        
        speczdiff = AtomPos_2_Cart[2]-AtomPos_1_Cart[2]
        print("LayerHeight:")
        print(speczdiff)
        return speczdiff
    
    def AtomDistance(self, TargetAtom_1, TargetAtom_2, p=True): # Atom2 - Atom1

        AtomPos_1_Cart, AtomPos_2_Cart = self._ReadPos(TargetAtom_1, TargetAtom_2, p=p)

        AtomDist_Cart = AtomPos_2_Cart - AtomPos_1_Cart
        specDist = np.sqrt(np.sum(AtomDist_Cart**2))
        print("Atom2 - Atom1:")
        print(AtomDist_Cart)
        # print("---")
        print("Distance:")
        print(specDist, "Å")
        return specDist

    def CrossSection_Rec(self, EndpointAtom_1, EndpointAtom_2, p=True):
        EndpointAtom_1_Cart, EndpointAtom_2_Cart = self._ReadPos(EndpointAtom_1, EndpointAtom_2, p=p)

        width  = np.sqrt(
                (EndpointAtom_2_Cart[0] - EndpointAtom_1_Cart[0])**2+(EndpointAtom_2_Cart[1]-EndpointAtom_1_Cart[1])**2
                )
        
        height = np.abs(EndpointAtom_2_Cart[2] - EndpointAtom_1_Cart[2])
        Area = width*height 

        print("Cross Section Area:")
        print(Area, "Å^2")

        return Area

    def CrossSection_Hex(self, EndpointAtom_1, EndpointAtom_2, p=True):
        EndpointAtom_1_Cart, EndpointAtom_2_Cart = self._ReadPos(EndpointAtom_1, EndpointAtom_2, p=p)
        length = (EndpointAtom_2_Cart[0] - EndpointAtom_1_Cart[0])**2+(EndpointAtom_2_Cart[1]-EndpointAtom_1_Cart[1])**2
        length = np.sqrt(length)
        
        Area = (3*np.sqrt(3)*(length**2))/2
        print("Cross Section Area:")
        print(Area, "Å^2")

        return Area
        

def File_Revise(JOBtype, R_Film, W_Film, Par_chang, Val_want):
    FileData = ""
    if JOBtype == 'Mumax3':
        with open(R_Film, 'r') as ReadF:
            for line in ReadF:
                if (Par_chang + ' := ') in line:
                    # line = line.replace(OldStr, NewStr)
                    line = Par_chang + ' := ' + str(Val_want) + '\n'
                elif (Par_chang + ' = ') in line:
                    line = Par_chang + ' = ' + str(Val_want) + '\n'
                FileData += line
    else:
        with open(R_Film, 'r') as ReadF:
            for line in ReadF:
                if Par_chang in line:
                    line = Val_want + '\n'

    with open(W_Film, 'w') as WriteF:
        WriteF.write(FileData)

def ColorList(p=False):
    color = ['#00107f', '#7f7f00', '#ff7a00', '#a000c8', '#ff0000', '#545659', '#79317b', 'green', '#373737']
    if p:
        plt.figure(figsize=(10, 1))
        for i, c in enumerate(color):
            plt.scatter(i, 0, color=c)
        plt.grid('--')
        plt.yticks([])
        plt.show()
    else:
        print("Loading color ...")

    return color

def Sum_Combination(n, Ntot):
    combinations = []
    for combo in combinations_with_replacement(range(1, Ntot+1), n):
        if sum(combo) == Ntot:
            combo = list(combo)
            combo.sort(reverse=True)
            combo = tuple(combo)
            combinations.append(combo)
    return combinations
