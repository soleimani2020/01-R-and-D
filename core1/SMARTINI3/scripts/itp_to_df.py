import numpy as np
import sys
import typing
import pandas as pd
from colors_text import TextColor as bcolors
import pandas as pd
import numpy as np
import sys
import typing
import os
import subprocess
from subprocess import call
import shutil
import multiprocessing
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from itertools import islice
import math



class ReadGro:
    """reading GRO file based on the doc"""

    info_msg: str = 'Message:\n'  # Message to pass for logging and writing
    line_len: int = 45  # Length of the lines in the data file
    gro_data: pd.DataFrame  # All the informations in the file
    # The follwings will set in __process_header_tail method:
    title: str  # Name of the system
    number_atoms: int  # Total number of atoms in the system
    pbc_box: str  # Size of the box (its 3 floats but save as a string)

    def __init__(self,
                 fname: str,  # Name of the input file
                 ) -> None:
        self.gro_data = self.read_gro(fname)

    def read_gro(self,
                 fname: str  # gro file name
                 ) -> pd.DataFrame:
        """read gro file lien by line"""
        counter: int = 0  # To count number of lines
        processed_line: list[dict[str, typing.Any]] = []  # All proccesed lines
        with open(fname, 'r', encoding='utf8') as f_r:
            while True:
                line = f_r.readline()
                if len(line) != self.line_len:
                    self.__process_header_tail(line.strip(), counter)
                else:
                    processed_line.append(self.__process_line(line.rstrip()))
                counter += 1
                if not line.strip():
                    break
        ReadGro.info_msg += f'\tFile name is {fname}\n'        
        ReadGro.info_msg += f'\tSystem title is {self.title}\n'
        ReadGro.info_msg += f'\tNumber of atoms is {self.number_atoms}\n'
        ReadGro.info_msg += f'\tBox boundary is {self.pbc_box}\n'
        return pd.DataFrame(processed_line)

    @staticmethod
    def __process_line(line: str  # Data line
                       ) -> dict[str, typing.Any]:
        """process lines of information"""
        resnr = int(line[0:5])
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        atomnr = int(line[15:20])
        a_x = float(line[20:28])
        a_y = float(line[28:36])
        a_z = float(line[36:44])
        processed_line: dict[str, typing.Any] = {
                                                 'residue_number': resnr,
                                                 'residue_name': resname,
                                                 'atom_name': atomname,
                                                 'atom_id': atomnr,
                                                 'x': a_x,
                                                 'y': a_y,
                                                 'z': a_z,
                                                }
        return processed_line

    def __process_header_tail(self,
                              line: str,  # Line in header or tail
                              counter: int  # Line number
                              ) -> None:
        """Get the header, number of atoms, and box size"""
        if counter == 0:
            self.title = line
        elif counter == 1:
            self.number_atoms = int(line)
        elif counter == self.number_atoms + 2:
            self.pbc_box = line




class APL_ANALYSIS:
    
    
    def __init__(self, membrane_LX: float = 43.51812, membrane_LY: float = 45.69877,  mesh_resolution: int = 1):  
        self.membrane_LX=membrane_LX
        self.membrane_LY=membrane_LY
        self.membrane_area = self.membrane_LX*self.membrane_LY
        self.mesh_resolution = mesh_resolution
        

    def _get_xy_grid(self) -> tuple[np.ndarray, np.ndarray]:
        """Generate a mesh grid for a given membrane area."""
        mesh_size_X = self.membrane_LX / self.mesh_resolution
        mesh_size_Y = self.membrane_LY / self.mesh_resolution
        grid_area=mesh_size_X*mesh_size_Y
        
        
        x_mesh, y_mesh = np.meshgrid(
            np.arange(0.0, self.membrane_LX, mesh_size_X),
            np.arange(0.0, self.membrane_LY, mesh_size_Y)
        )

        Mesh_NUMBER=1
        

        return x_mesh, y_mesh , grid_area ,  Mesh_NUMBER , self.mesh_resolution , mesh_size_X , mesh_size_Y  , self.membrane_LX  , self.membrane_LY
    
    

    @staticmethod
    def write_gromacs_gro(gro_data: pd.DataFrame,
                          filename: str,  # Name of the output file
                          pbc_box=None,
                          title=None
                          ) -> None:
        """Write DataFrame to a GROMACS gro file."""
        
        
        df_i: pd.DataFrame = gro_data.copy()
        
        output_file_path = os.path.join(filename)
        
        with open(output_file_path, 'w', encoding='utf8') as gro_file:
            if title:
                gro_file.write(f'{title}')  # Add a comment line
            gro_file.write(f'{len(df_i)}\n')  # Write the total number of atoms
            for _, row in df_i.iterrows():
                line = f'{row["residue_number"]:>5}' \
                       f'{row["residue_name"]:<5}' \
                       f'{row["atom_name"]:>5}' \
                       f'{row["atom_id"]:>5}' \
                       f'{row["x"]:8.3f}' \
                       f'{row["y"]:8.3f}' \
                       f'{row["z"]:8.3f}\n'
                gro_file.write(line)
            if pbc_box:
                gro_file.write(f'{pbc_box}\n')


    
    @classmethod
    def process_mesh(cls, x_mesh, y_mesh, mesh_size_X, mesh_size_Y, Mesh_NUMBER, mesh_resolution, xyz_i, max_z_threshold, min_z_threshold, frame):
        selected_atoms_info = {}

        selected_atoms_data_original = xyz_i
        column_names_original = ['residue_number', 'residue_name', 'atom_name', 'atom_id', 'x', 'y', 'z']
        df_selected_Raws_original = pd.DataFrame(selected_atoms_data_original, columns=column_names_original)



    
        # Define the conditions for selecting atoms
        x_min_mesh = -100
        x_max_mesh = 1.518
        y_min_mesh = 1.3
        y_max_mesh = 100

        #mask_nc3 = np.full(len(xyz_i), True)
        #ind_in_mesh_nc3_before_mask = np.arange(len(xyz_i))

        #mask_nc3 = (xyz_i[:, 2] == 'BB')
    
        # Apply the mask to select atoms within the current mesh element based on XY & Z
        ind_in_mesh_nc3 = np.where(
            (xyz_i[:, 4] >= x_min_mesh) &
            (xyz_i[:, 4] < x_max_mesh) &
            (xyz_i[:, 5] >= y_min_mesh) &
            (xyz_i[:, 5] < y_max_mesh) &
            (xyz_i[:, 6] < max_z_threshold) &
            (xyz_i[:, 6] > min_z_threshold) )
    
        # Select atom data based on the indices obtained
        selected_atoms_data = xyz_i[ind_in_mesh_nc3]
        #print(selected_atoms_data)
    
        column_names = ['residue_number', 'residue_name', 'atom_name', 'atom_id', 'x', 'y', 'z']
        df_selected_Raws = pd.DataFrame(selected_atoms_data, columns=column_names)
        #print(len(df_selected_Raws))
        #print(df_selected_Raws)

        Index_List = df_selected_Raws.iloc[:, 3].tolist()




        
        # Function to calculate distance between two particles
        def calculate_distance(particle1, particle2):
            x1, y1, z1 = particle1['x'], particle1['y'], particle1['z']
            x2, y2, z2 = particle2['x'], particle2['y'], particle2['z']
            return math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
        
        # Iterate through all combinations of particles and calculate distances
        distances = []
        for i in range(len(df_selected_Raws)):
            for j in range(i + 1, len(df_selected_Raws)):
                particle1 = df_selected_Raws.iloc[i]
                particle2 = df_selected_Raws.iloc[j]
                
                # Example: Update particle numbers with a prefix 'A' and 'B'
                updated_particle1 = f'A{particle1["residue_number"]}'
                updated_particle2 = f'B{particle2["residue_number"]}'
                
                distance = calculate_distance(particle1, particle2)
                distances.append({
                    'particle1': updated_particle1,
                    'particle2': updated_particle2,
                    'distance': distance
                })
        
        # Create a DataFrame from the distances list
        df_distances_aftermask = pd.DataFrame(distances)
        #print("df_distances_aftermask:",df_distances_aftermask)
        


        # Filter rows for particles with indices in Index_List
        df_selected_Raws_filtered = df_selected_Raws_original[df_selected_Raws_original['atom_id'].isin(Index_List)]
        
        # Iterate through all combinations of particles and calculate distances
        distances = []
        for i in range(len(df_selected_Raws_filtered)):
            for j in range(i + 1, len(df_selected_Raws_filtered)):
                particle1 = df_selected_Raws_filtered.iloc[i]
                particle2 = df_selected_Raws_filtered.iloc[j]
                
                # Update particle numbers if needed
                updated_particle1 = particle1['residue_number']
                updated_particle2 = particle2['residue_number']
                
                distance = calculate_distance(particle1, particle2)
                distances.append({
                    'particle1': updated_particle1,
                    'particle2': updated_particle2,
                    'distance': distance
                })
        
        # Create a DataFrame from the distances list
        df_distances_beforemask = pd.DataFrame(distances)
        #print("df_distances_beforemask:",df_distances_beforemask)


        # Check if the last columns are exactly equal
        last_column_aftermask = df_distances_aftermask.iloc[:, -1]
        last_column_beforemask = df_distances_beforemask.iloc[:, -1]
        
        #if last_column_aftermask.equals(last_column_beforemask):
            #print("Distances are exactly equal.")
        #else:
            #print("Something wrong ! Distances are not exactly equal.")
                
        

        
        apl_instance = APL_ANALYSIS()
        filename = "NEW_DYNAMINE.gro"  
        pbc_box = "43.90298 46.10292 17.86046"
        title = "This file contains head part of the dynamine.\n"
        
        # Write the DataFrame to a GROMACS GRO file without specifying the output directory
        apl_instance.write_gromacs_gro(df_selected_Raws, filename, pbc_box, title)

        #print("Indexes in gro section\n:",Index_List)
    
        return Index_List
        



##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################


# A helper function needed by most of the classes to clean the lines
def free_char_line(line: str  # line of the itp file
                   ) -> list[str]:  # Free from chars
    """cheack the lines and return the line free special chars"""
    char_list: list[str] = [';', '#', ':', '...']  # chars eliminate in lines
    l_line: list[str]  # Breaking the line cahrs
    l_line = line.strip().split(' ')
    l_line = [item for item in l_line if item]
    l_line = [item for item in l_line if item not in char_list]
    return l_line


# A helper function needed by most of the classes to get types for LAMMPS
def get_type(lst: list[str]  # list to get the number of distenguished ones
             ) -> list[int]:  # types' index
    """make type based on the unique items in the lst"""
    type_set: set[str] = set(lst)  # eleminate the repeated names
    type_dict: dict[str, int]  # to make a list with type index
    type_dict = {item: i+1 for i, item in enumerate(type_set)}
    types: list[int]  # types to return
    types = [type_dict[item] for item in lst]
    return types


class Itp:
    """read itp file (protein.itp) and return a DataFrame of the information
    within the file"""
    def __init__(self,
                 fname: str  # Name of the itp file
                 ) -> None:
        self.get_itp(fname)

    def get_itp(self,
                fname: str  # Name of the itp file
                ) -> None:
        """read the file line by line and call related methods"""
        atoms: bool = False  # Flag of 'atoms' occurrence
        bonds: bool = False  # Flag of 'bonds' occurrence
        constraints: bool = False  # Flag of 'constraints' occurrence
        angles: bool = False  # Flag of 'angles' occurrence
        dihedrals: bool = False  # Flag of 'dihedrals' occurrence
        impropers: bool = False  # Flag of 'impropers' occurrence
        moleculetype: bool = False  # Flag of 'moleculetype' occurrence
        atoms_info: list[str] = []  # to append atoms lines
        bonds_info: list[str] = []  # to append bonds lines
        angles_info: list[str] = []  # to append angles lines
        dihedrals_info: list[str] = []  # to append dihedrals lines

        constraints_info : list[str] = []  # to append dihedrals lines
        
        with open(fname, 'r') as f:
            while True:
                line: str = f.readline()
                if line.strip():
                    if line.strip().startswith('['):
                        wilds: list[str]  # Parts of the line
                        wilds = line.strip().split()
                        if wilds[1] == 'atoms':
                            atoms, bonds, angles, dihedrals, impropers,\
                                moleculetype, constraints = True, False, False, False,\
                                False, False, False
                        elif wilds[1] == 'bonds':
                            atoms, bonds, angles, dihedrals, impropers,\
                                moleculetype, constraints = False, True, False, False,\
                                False, False, False
                        elif wilds[1] == 'angles':
                            atoms, bonds, angles, dihedrals, impropers,\
                                moleculetype, constraints = False, False, True, False,\
                                False, False, False
                        elif wilds[1] == 'dihedrals':
                            atoms, bonds, angles, dihedrals, impropers,\
                                moleculetype , constraints= False, False, False, True,\
                                False, False, False
                        elif wilds[1] == 'moleculestype':
                            atoms, bonds, angles, dihedrals, impropers,\
                                moleculetype, constraints = False, False, False, False,\
                                False, True, False
                        elif wilds[1] == 'constraints':
                            atoms, bonds, angles, dihedrals, impropers,\
                                moleculetype , constraints  = False, False, False, False,\
                                False, False , True 

                        else:
                            atoms, bonds, angles, dihedrals, impropers,\
                                moleculetype , constraints  = False, False, False, False,\
                                False, False , False 
                    
                    else:
                        if atoms:
                            atoms_info.append(line)
                        if bonds:
                            bonds_info.append(line)
                        if angles:
                            angles_info.append(line)
                        if dihedrals:
                            if 'impropers' not in free_char_line(line):
                                impropers = False
                            else:
                                impropers = True
                                dihedrals = False
                        if dihedrals and not impropers:
                            dihedrals_info.append(line)
                        if constraints :
                            constraints_info.append(line)
                            

                
                if not line:
                    break
                


        #print(constraints_info)
        atom = AtomsInfo(atoms_info)

        self.exchange_dict=atom.exchange_dict
        
        bond = BondsInfo(atoms=atom.df, bonds=bonds_info,exchange_dict=self.exchange_dict)
        angle = AnglesInfo(atoms=atom.df, angles=angles_info,exchange_dict=self.exchange_dict)
        dihedral = DihedralsInfo(atoms=atom.df, dihedrals=dihedrals_info,exchange_dict=self.exchange_dict)
        constraint= ConstraintsInfo(atoms=atom.df, constraints=constraints_info,exchange_dict=self.exchange_dict)
        
        self.atoms_extra = atom.df
        self.bonds = bond.df
        self.angles = angle.df
        self.dihedrals = dihedral.df
        self.constraints = constraint.df
        

        

        


####### AtomsInfo  class ##################################################################################################


class AtomsInfo:
    """get atoms wild information and retrun a DataFrame"""
    def __init__(self,
                 atoms: list[str]  # lines read by Itp class
                 ) -> None:
        self.df , self.exchange_dict = self.get_atoms_info(atoms)

    def get_atoms_info(self,
                       atoms: list[str]  # Lines of the atoms' section
                       ) -> pd.DataFrame:
        """get atoms info from the file"""
        l_line: list[str]  # Breaking the line cahrs
        # Check if header of the atoms section is same as the defeined one
        columns: list[str]   # columns for the atoms dict, name of each column
        columns_al: list[str]   # Alternative columns names
        columns = ['atomnr', 'atomtype', 'resnr', 'resname', 'atomname',
                   'chargegrp', 'charge', 'mass']
        columns_al = ['nr', 'type', 'resnr', 'resid', 'atom',
                    'cgnr', 'charge', 'mass']
        atomnr: list[typing.Any] = []  # list to append info: atoms id
        atomtype: list[typing.Any] = []  # list to append info: forcefield type
        resnr: list[typing.Any] = []  # list to append info: res infos
        resname: list[typing.Any] = []  # list to append info: res number
        atomname: list[typing.Any] = []  # list to append info: atom name
        chargegrp: list[typing.Any] = []  # list to append info: charge group
        charge: list[typing.Any] = []  # list to append info: charge value
        mass: list[typing.Any] = []  # list to append info: mass value
        # The 3 columns below do not have header, I added to get useful name
        atomsty: list[typing.Any] = []  # list to append info: name with style
        chemi: list[typing.Any] = []  # list to append info: name with alkane
        name: list[typing.Any] = []  # list to append info: real name
        columns_falg: bool = False  # if 10 columns data
        columns_al_flag: bool = False  # if 8 columns data
        for line in atoms:
            if line.startswith(';'):  # line start with ';' are commets&header
                l_line = free_char_line(line)
                if 'Total' not in l_line:  # Not header!
                    if l_line == columns:
                        columns_falg = True
                    elif l_line == columns_al:
                        columns_al_flag = True
                    else:
                        exit(f'{bcolors.FAIL}{self.__class__.__name__}:\n'
                             f'\tError in the [ atoms ] header of the '
                             f'itp file\n{bcolors.ENDC}')
            else:
                l_line = free_char_line(line)
                atomnr.append(l_line[0])
                atomtype.append(l_line[1])
                resnr.append(l_line[2])
                resname.append(l_line[3])
                atomname.append(l_line[4])
                chargegrp.append(l_line[5])
                charge.append(l_line[6])
                mass.append(l_line[7])
                if columns_falg:
                    atomsty.append(l_line[8])
                    chemi.append(l_line[9])
                    name.append(l_line[10])
        df: pd.DataFrame  # DataFrame from the infos
        df = pd.DataFrame(columns=columns)
        
        df['atomnr'] = atomnr
        df['atomtype'] = atomtype
        df['resnr'] = resnr
        df['resname'] = resname
        df['atomname'] = atomname
        df['chargegrp'] = chargegrp
        df['charge'] = charge
        df['mass'] = mass

        #print(df)


        # Drop rows based on indices

        frame=0
        read_gro_instance = ReadGro(fname="confout.gro")
        gro_data = read_gro_instance.gro_data
        xyz_i = gro_data[['residue_number', 'residue_name', 'atom_name', 'atom_id','x', 'y', 'z']].values
        mesh_generator = APL_ANALYSIS()
        x_mesh, y_mesh , grid_area ,  Mesh_NUMBER , mesh_resolution , mesh_size_X , mesh_size_Y , membrane_LX , membrane_LY= mesh_generator._get_xy_grid()
        Indices = mesh_generator.process_mesh(x_mesh, y_mesh, mesh_size_X, mesh_size_Y, Mesh_NUMBER, mesh_resolution,  xyz_i=xyz_i, max_z_threshold=1000, min_z_threshold=-1000, frame=frame)

        #print("Indices in atom section:\n", Indices)

        Indices = [index - 1 for index in Indices]

        
        df = df.loc[Indices]

        #print(df)

        atomnr_list_before = df['atomnr'].tolist()
        df['atomnr'] = range(1, len(df) + 1)
        #df['resnr'] = range(1, len(df) + 1)
        df['chargegrp'] = range(1, len(df) + 1)
        value_mapping = {value: index for index, value in enumerate(range(1, len(df)+1), start=1)}
        df['atomnr'] = df['atomnr'].astype(int)
        df['atomnr'] = df['atomnr'].map(value_mapping)
        #df['resnr'] = df['resnr'].astype(int)
        #df['resnr'] = df['resnr'].map(value_mapping)
        df['chargegrp'] = df['chargegrp'].astype(int)
        df['chargegrp'] = df['chargegrp'].map(value_mapping)
        atomnr_list_after = df['atomnr'].tolist()
        exchange_dict = dict(zip(atomnr_list_before, atomnr_list_after))
        #print("exchange_dict in atom section:\n",exchange_dict)

        #print(df)

        
        # Create a NumPy array
        ATOMS = np.zeros([len(df), 7], dtype=object)
        for i in range(7):
            for j in range(len(df)):
                ATOMS[j][i] = df.iloc[j, i]  # Assuming i is the correct column index
        
        # Save the NumPy array to a file
        np.savetxt("New_Atom.itp", ATOMS, delimiter=" ", fmt="%8s %6s %4s %7s %7s %7s %7s  ")
        
    
        return df , exchange_dict






####### BondsInfo  class ##################################################################################################


class BondsInfo:



    """get the bonds list from Itp class and return a dataframe"""
    def __init__(self,
                 bonds: list[str],  # lines of bonds section read by Itp class
                 atoms: pd.DataFrame,  # atoms df from AtomsInfo to get names
                 exchange_dict
                 ) -> None:
        """get the bonds infos"""
        self.exchange_dict=exchange_dict
        
        self.mk_bonds_df(bonds, atoms)
        

    def mk_bonds_df(self,
                    bonds: list[str],  # lines of bonds section
                    atoms: pd.DataFrame , # atoms df from AtomInfo
                    ) -> None:
        """call all the methods to make the bonds DataFrame"""
        ai: list[int]  # index of the 1st atoms in the bonds
        aj: list[int]  # index of the 2nd atoms in the bonds
        names: list[str]  # name of the bonds
        ai, aj, typ , r0 , k0 , names = self.get_bonds(bonds)
        self.df = self.mk_df( ai, aj, typ , r0 , k0 , names, atoms)

    def get_bonds(self,
                  bonds: list[str]  ,   # lines of bonds section read by Itp class
                  ) -> pd.DataFrame:  # DataFrame of the bonds
        """return bonds dataframe to make bonds dataframe"""
        columns: list[str]  # Columns of the bonds wild
        columns = ['ai', 'aj', 'funct', 'r', 'k']
        columns_al =  ['ai', 'aj', 'fu']
        ai: list[int] = []  
        aj: list[int] = [] 
        typ: list[int] = []  
        r0: list[int] = []  
        k0: list[int] = []  
        
        names: list[str] = []  # name of the bonds
        for line in bonds:
            if line.startswith(';'):  # line start with ';' are commets&header
                l_line = free_char_line(line)
                if 'Total' not in l_line:  # Not header!
                    if l_line == columns or l_line == columns_al:
                        pass
                    else:
                        exit(f'{bcolors.FAIL}{self.__class__.__name__}:\n'
                             f'\tError in the [ bonds ] header of the '
                             f'itp file\n{bcolors.ENDC}')
            else:
                l_line = free_char_line(line)
                ai.append(int(l_line[0]))
                aj.append(int(l_line[1]))
                typ.append(int(l_line[2]))
                r0.append(float(l_line[3]))
                k0.append(float(l_line[4]))
                names.append("l_line[4]")
        return ai, aj, typ , r0 , k0 , names

    def mk_df(self,
              ai, aj, typ , r0 , k0 , names , 
              atoms: pd.DataFrame  # atoms df from AtomsInfo to cehck the name
              ) -> pd.DataFrame:   # bonds DataFrame


        """make DataFrame"""
        
        df: pd.DataFrame  
        df = pd.DataFrame(columns=['ai', 'aj', 'typ', 'r0', 'k0'])
        df['ai'] = ai
        df['aj'] = aj
        df['r0'] = r0
        df['k0'] = k0
        df['typ'] = typ
        df.index += 1



        # Drop rows based on indices
        frame=0
        read_gro_instance = ReadGro(fname="confout.gro")
        gro_data = read_gro_instance.gro_data
        xyz_i = gro_data[['residue_number', 'residue_name', 'atom_name', 'atom_id','x', 'y', 'z']].values
        mesh_generator = APL_ANALYSIS()
        x_mesh, y_mesh , grid_area ,  Mesh_NUMBER , mesh_resolution , mesh_size_X , mesh_size_Y , membrane_LX , membrane_LY= mesh_generator._get_xy_grid()
        Indices = mesh_generator.process_mesh(x_mesh, y_mesh, mesh_size_X, mesh_size_Y, Mesh_NUMBER, mesh_resolution,  xyz_i=xyz_i, max_z_threshold=1000, min_z_threshold=-1000, frame=frame)
        #print("Indice for bonds:\n",Indices) 
        #Indices = [index + 1 for index in Indices]
        #print("Indice for bonds:\n",Indices) 

        # Keep rows where 'ai' and 'aj' is in the list Indices
        mask = (df['ai'].isin(Indices)) & (df['aj'].isin(Indices))
        df = df[mask]


        Graphene = np.zeros([len(df), 5], dtype=object)
        
        for i in range(5):
            for j in range(len(df)):
                # Use the correct DataFrame column names or indices
                Graphene[j][i] = df.iloc[j, i]  # Assuming i is the correct column index
        
        # Convert the first three columns to integers
        Graphene[:, :3] = Graphene[:, :3].astype(int)
        
        # Save the NumPy array to a file
        np.savetxt("OLD_bond.itp", Graphene, delimiter=" ", fmt="%8d %6d %4d %7s %7s")
                
        #print(df)

        # Convert the 'ai' and 'aj'  columns to string type if they're not already
        df['ai'] = df['ai'].astype(str)
        df['aj'] = df['aj'].astype(str)
        
        # Update 'ai' and 'aj' columns using the exchange_dict
        df['ai'] = df['ai'].map(self.exchange_dict).fillna(df['ai'])
        df['aj'] = df['aj'].map(self.exchange_dict).fillna(df['aj'])

        #print(df)


        # Create a NumPy array
        Graphene = np.zeros([len(df), 5], dtype=object)
        
        for i in range(5):
            for j in range(len(df)):
                # Use the correct DataFrame column names or indices
                Graphene[j][i] = df.iloc[j, i]  # Assuming i is the correct column index
        
        # Convert the first three columns to integers
        Graphene[:, :3] = Graphene[:, :3].astype(int)
        
        # Save the NumPy array to a file
        np.savetxt("New_bond.itp", Graphene, delimiter=" ", fmt="%8d %6d %4d %7s %7s")
        
        return df



####### ConstraintsInfo  class ##################################################################################################


class ConstraintsInfo:



    """get the consts list from Itp class and return a dataframe"""
    def __init__(self,
                 constraints: list[str],  # lines of consts section read by Itp class
                 atoms: pd.DataFrame,  # atoms df from AtomsInfo to get names
                 exchange_dict
                 ) -> None:
        """get the consts infos"""
        self.exchange_dict=exchange_dict
        
        self.mk_consts_df(constraints, atoms)
        

    def mk_consts_df(self,
                    consts: list[str],  # lines of consts section
                    atoms: pd.DataFrame , # atoms df from AtomInfo
                    ) -> None:
        """call all the methods to make the consts DataFrame"""
        ai: list[int]  # index of the 1st atoms in the consts
        aj: list[int]  # index of the 2nd atoms in the consts
        names: list[str]  # name of the consts
        ai, aj, typ , r0  , names = self.get_consts(consts)
        self.df = self.mk_df( ai, aj, typ , r0  , names, atoms)

    def get_consts(self,
                  consts: list[str]  ,   # lines of consts section read by Itp class
                  ) -> pd.DataFrame:  # DataFrame of the consts
        """return consts dataframe to make consts dataframe"""
        columns: list[str]  # Columns of the consts wild
        columns = ['ai', 'aj', 'funct', 'r']
        columns_al =  ['ai', 'aj', 'fu']
        ai: list[int] = []  
        aj: list[int] = [] 
        typ: list[int] = []  
        r0: list[int] = []  


        
        names: list[str] = []  # name of the bonds
        for line in consts:
            if line.startswith(';'):  # line start with ';' are commets&header
                l_line = free_char_line(line)
                if 'Total' not in l_line:  # Not header!
                    if l_line == columns or l_line == columns_al:
                        pass
                    else:
                        exit(f'{bcolors.FAIL}{self.__class__.__name__}:\n'
                             f'\tError in the [ bonds ] header of the '
                             f'itp file\n{bcolors.ENDC}')
            else:
                l_line = free_char_line(line)
                ai.append(int(l_line[0]))
                aj.append(int(l_line[1]))
                typ.append(int(l_line[2]))
                r0.append(float(l_line[3]))
                names.append("consts")


        print("Constraint Section:")
        #print(ai, aj, typ , r0  , names, end="\n")
        return ai, aj, typ , r0  , names

    def mk_df(self,
              ai, aj, typ , r0  , names , 
              atoms: pd.DataFrame  # atoms df from AtomsInfo to cehck the name
              ) -> pd.DataFrame:   # consts DataFrame


        """make DataFrame"""
        
        df: pd.DataFrame  
        df = pd.DataFrame(columns=['ai', 'aj', 'typ', 'r0'])
        df['ai'] = ai
        df['aj'] = aj
        df['r0'] = r0
        df['typ'] = typ
        df.index += 1

        df.to_csv("constraints",sep=" ")



        # Drop rows based on indices
        frame=0
        read_gro_instance = ReadGro(fname="confout.gro")
        gro_data = read_gro_instance.gro_data
        xyz_i = gro_data[['residue_number', 'residue_name', 'atom_name', 'atom_id','x', 'y', 'z']].values
        mesh_generator = APL_ANALYSIS()
        x_mesh, y_mesh , grid_area ,  Mesh_NUMBER , mesh_resolution , mesh_size_X , mesh_size_Y , membrane_LX , membrane_LY= mesh_generator._get_xy_grid()
        Indices = mesh_generator.process_mesh(x_mesh, y_mesh, mesh_size_X, mesh_size_Y, Mesh_NUMBER, mesh_resolution,  xyz_i=xyz_i, max_z_threshold=1000, min_z_threshold=-1000, frame=frame)
        print("Indice for constraints:\n",Indices) 

        

        # Keep rows where 'ai' and 'aj' is in the list Indices
        mask = (df['ai'].isin(Indices)) & (df['aj'].isin(Indices))
        df = df[mask]


        Graphene = np.zeros([len(df), 4], dtype=object)
        
        for i in range(4):
            for j in range(len(df)):
                # Use the correct DataFrame column names or indices
                Graphene[j][i] = df.iloc[j, i]  # Assuming i is the correct column index
        
        # Convert the first three columns to integers
        Graphene[:, :2] = Graphene[:, :2].astype(int)
        
        # Save the NumPy array to a file
        np.savetxt("OLD_consts.itp", Graphene, delimiter=" ", fmt="%8d %6d %4d %7s")
                
        print(df)

        # Convert the 'ai' and 'aj'  columns to string type if they're not already
        df['ai'] = df['ai'].astype(str)
        df['aj'] = df['aj'].astype(str)
        
        # Update 'ai' and 'aj' columns using the exchange_dict
        df['ai'] = df['ai'].map(self.exchange_dict).fillna(df['ai'])
        df['aj'] = df['aj'].map(self.exchange_dict).fillna(df['aj'])

        print(df)



        # Create a NumPy array
        Graphene = np.zeros([len(df), 4], dtype=object)
        
        for i in range(4):
            for j in range(len(df)):
                # Use the correct DataFrame column names or indices
                Graphene[j][i] = df.iloc[j, i]  # Assuming i is the correct column index
        
        # Convert the first four columns to integers
        Graphene[:, :3] = Graphene[:, :3].astype(int)
        
        # Save the NumPy array to a file
        np.savetxt("New_constraint.itp", Graphene, delimiter=" ", fmt="%8d %6d %4d %7s")



        
        return df




####### AnglesInfo  class ##################################################################################################


class AnglesInfo:
    """get the angles list from Itp class and return a dataframe"""
    def __init__(self,
                 angles: list[str],  # lines of angles section by Itp class
                 atoms: pd.DataFrame,  # atoms df from AtomsInfo to get names
                 exchange_dict
                 ) -> None:
        """get the angles infos"""
        self.exchange_dict=exchange_dict
        
        self.mk_angles_df(angles, atoms)

    def mk_angles_df(self,
                     angles: list[str],  # lines of angles section
                     atoms: pd.DataFrame  # atoms df from AtomInfo
                     ) -> None:
        """call all the methods to make the bonds DataFrame"""
        ai: list[int]  # index of the 1st atoms in the angles
        aj: list[int]  # index of the 2nd atoms in the angles
        ak: list[int]  # index of the 3rd atoms in the angles
        names: list[str]  # name of the angles
        ai, aj, ak, typ, r0 , k0 , names = self.get_angles(angles)
        self.df = self.mk_df(ai, aj, ak, typ , r0 , k0 , names, atoms)

    def get_angles(self,
                   angles: list[str],  # lines of angles section by Itp class
                   ) -> pd.DataFrame:  # DataFrame of the angles
        """return bonds dataframe to make angles dataframe"""
        columns: list[str]  # Columns of the angles wild
        columns = ['ai', 'aj', 'ak', 'funct', 'theta', 'cth']
        ai: list[int] = []  # index of the 1st atoms in the angles
        aj: list[int] = []  # index of the 2nd atoms in the angles
        ak: list[int] = []  # index of the 3rd atoms in the angles
        typ: list[int] = []  # index of the 2nd atoms in the bonds
        r0: list[int] = []  # index of the 2nd atoms in the bonds
        k0: list[int] = []  # index of the 2nd atoms in the bonds


        
        names: list[str] = []  # name of the angles
        columns_al =  ['ai', 'aj', 'ak', 'fu']

        for line in angles:
            if line.startswith(';'):  # line start with ';' are commets&header
                l_line = free_char_line(line)
                if 'Total' not in l_line:  # Not header!
                    if l_line == columns:
                        pass
                    else:
                        exit(f'{bcolors.FAIL}{self.__class__.__name__}:\n'
                             f'\tError in the [ angles ] header of the '
                             f'itp file\n{bcolors.ENDC}')
            else:
                l_line = free_char_line(line)
                ai.append(int(l_line[0]))
                aj.append(int(l_line[1]))
                ak.append(int(l_line[2]))
                typ.append(int(l_line[3]))
                r0.append(float(l_line[4]))
                k0.append(float(l_line[5]))
                names.append("l_line[6]")
        return ai, aj, ak, typ , r0 , k0 , names

    def mk_df(self,
              ai, aj, ak, typ , r0 , k0 , names , 
              atoms: pd.DataFrame  # atoms df from AtomsInfo to cehck the name
              ) -> pd.DataFrame:  # angles DataFrame
        """make DataFrame and check if they are same as atoms name"""
        df: pd.DataFrame  # to save the angles_df
        df = pd.DataFrame(columns=['ai', 'aj', 'ak','typ', 'r0', 'k0'])
        df['ai'] = ai
        df['aj'] = aj
        df['ak'] = ak
        df['r0'] = r0
        df['k0'] = k0
        df['typ'] = typ
        df.index += 1
        # Drop rows based on indices

        # Drop rows based on indices

        frame=0
        read_gro_instance = ReadGro(fname="confout.gro")
        gro_data = read_gro_instance.gro_data
        xyz_i = gro_data[['residue_number', 'residue_name', 'atom_name', 'atom_id','x', 'y', 'z']].values
        mesh_generator = APL_ANALYSIS()
        x_mesh, y_mesh , grid_area ,  Mesh_NUMBER , mesh_resolution , mesh_size_X , mesh_size_Y , membrane_LX , membrane_LY= mesh_generator._get_xy_grid()
        Indices = mesh_generator.process_mesh(x_mesh, y_mesh, mesh_size_X, mesh_size_Y, Mesh_NUMBER, mesh_resolution,  xyz_i=xyz_i, max_z_threshold=1000, min_z_threshold=-1000, frame=frame)
        
        #print("Indices in the angle section\n:",Indices)

        
        # Keep rows where 'ai' and  'aj' and 'ak' is in the list Indices
        mask = (df['ai'].isin(Indices)) & (df['aj'].isin(Indices)) & (df['ak'].isin(Indices))
        df = df[mask]


        Graphene = np.zeros([len(df), 6], dtype=object)
        
        for i in range(6):
            for j in range(len(df)):
                # Use the correct DataFrame column names or indices
                Graphene[j][i] = df.iloc[j, i]  # Assuming i is the correct column index
        
        # Convert the first four columns to integers
        Graphene[:, :4] = Graphene[:, :4].astype(int)
        
        # Save the NumPy array to a file
        np.savetxt("OLD_angle.itp", Graphene, delimiter=" ", fmt="%8d %6d %4d %7s %7s %7s")
        

        #print(df)

        
        # Convert the 'ai' and 'aj' columns to string type if they're not already
        df['ai'] = df['ai'].astype(str)
        df['aj'] = df['aj'].astype(str)
        df['ak'] = df['ak'].astype(str)
        
        # Update 'ai' and 'aj', 'ak' columns using the exchange_dict
        df['ai'] = df['ai'].map(self.exchange_dict).fillna(df['ai'])
        df['aj'] = df['aj'].map(self.exchange_dict).fillna(df['aj'])
        df['ak'] = df['ak'].map(self.exchange_dict).fillna(df['ak'])

        #print(df)



        # Create a NumPy array
        Graphene = np.zeros([len(df), 6], dtype=object)
        
        for i in range(6):
            for j in range(len(df)):
                # Use the correct DataFrame column names or indices
                Graphene[j][i] = df.iloc[j, i]  # Assuming i is the correct column index
        
        # Convert the first four columns to integers
        Graphene[:, :4] = Graphene[:, :4].astype(int)
        
        # Save the NumPy array to a file
        np.savetxt("New_angle.itp", Graphene, delimiter=" ", fmt="%8d %6d %4d %7s %7s %7s")
        return df


####### DihedralsInfo  class ##################################################################################################



class DihedralsInfo:
    """get the dihdrals list from Itp class and return a dataframe"""
    def __init__(self,
                 dihedrals: list[str],  # lines of dihedrals section by Itp
                 atoms: pd.DataFrame , # atoms df from AtomsInfo to get names
                 exchange_dict
                 ) -> None:
        """get the dihedrals infos"""
        self.exchange_dict=exchange_dict
        self.mk_dihedrals_df(dihedrals, atoms)

    def mk_dihedrals_df(self,
                        dihedrals: list[str],  # lines of dihedrals section
                        atoms: pd.DataFrame  # atoms df from AtomInfo
                        ) -> None:
        """call all the methods to make the bonds DataFrame"""
        ai: list[int]  # index of the 1st atoms in the dihedrals
        aj: list[int]  # index of the 2nd atoms in the dihedrals
        ak: list[int]  # index of the 3rd atoms in the dihedrals
        ah: list[int]  # index of the 4th atoms in the dihedrals
        names: list[str]  # name of the dihedrals
        ai, aj, ak, al, typ , r0 , k0 , names = self.get_dihedrals(dihedrals)
        self.df = self.mk_df(ai, aj, ak, al, typ , r0 , k0 , names, atoms)

    def get_dihedrals(self,
                      dihedrals: list[str],  # lines of dihedrals section
                      ) -> pd.DataFrame:  # DataFrame of the dihedrals
        """return bonds dataframe to make dihedrals dataframe"""
        columns: list[str]  # Columns of the dihedrals wild
        columns = ['ai', 'aj', 'ak', 'al', 'funct', 'C0', 'C5']
        ai: list[int] = []  # index of the 1st atoms in the dihedrals
        aj: list[int] = []  # index of the 2nd atoms in the dihedrals
        ak: list[int] = []  # index of the 3rd atoms in the dihedrals
        ah: list[int] = []  # index of the 4th atoms in the dihedrals
        typ: list[int] = []  # index of the 2nd atoms in the bonds
        r0: list[int] = []  # index of the 2nd atoms in the bonds
        k0: list[int] = []  # index of the 2nd atoms in the bonds
        names: list[str] = []  # name of the dihedrals
        for line in dihedrals:
            if line.startswith(';'):  # line start with ';' are commets&header
                l_line = free_char_line(line)
                if 'Total' not in l_line:  # Not header!
                    if l_line == columns:
                        pass
                    else:
                        exit(f'{bcolors.FAIL}{self.__class__.__name__}:\n'
                             f'\tError in the [ dihedrals ] header of the '
                             f'itp file\n{bcolors.ENDC}')
            else:
                l_line = free_char_line(line)
                ai.append(int(l_line[0]))
                aj.append(int(l_line[1]))
                ak.append(int(l_line[2]))
                ah.append(int(l_line[3]))
                typ.append(int(l_line[4]))
                r0.append(float(l_line[5]))
                k0.append(float(l_line[6]))
                names.append("l_line[7]")

        return ai, aj, ak, ah, typ , r0 , k0 , names

    def mk_df(self,
              ai, aj, ak, ah, typ , r0 , k0 , names , 
              atoms: pd.DataFrame  # atoms df from AtomsInfo to cehck the name
              ) -> pd.DataFrame:  # dihedrals DataFrame
        """make DataFrame and check if they are same as atoms name"""
        df: pd.DataFrame  # to save the dihedrals_df
        df = pd.DataFrame(
            columns=[ 'ai', 'aj', 'ak', 'ah', 'typ','r0', 'k0'])
        df['ai'] = ai
        df['aj'] = aj
        df['ak'] = ak
        df['ah'] = ah
        df['r0'] = r0
        df['k0'] = k0
        df['typ'] = typ
        df.index += 1


        
        # Drop rows based on indices


        # Drop rows based on indices

        frame=0
        read_gro_instance = ReadGro(fname="confout.gro")
        gro_data = read_gro_instance.gro_data
        xyz_i = gro_data[['residue_number', 'residue_name', 'atom_name', 'atom_id','x', 'y', 'z']].values
        mesh_generator = APL_ANALYSIS()
        x_mesh, y_mesh , grid_area ,  Mesh_NUMBER , mesh_resolution , mesh_size_X , mesh_size_Y , membrane_LX , membrane_LY= mesh_generator._get_xy_grid()
        Indices = mesh_generator.process_mesh(x_mesh, y_mesh, mesh_size_X, mesh_size_Y, Mesh_NUMBER, mesh_resolution,  xyz_i=xyz_i, max_z_threshold=1000, min_z_threshold=-1000, frame=frame)

        #print("Indices in the dihedral section:\n",Indices)

        
        # Keep rows where 'ai' or 'aj' or 'ak' or 'ah' is in the list Indices
        mask = (df['ai'].isin(Indices)) & (df['aj'].isin(Indices)) & (df['ak'].isin(Indices)) & (df['ah'].isin(Indices))
        df = df[mask]



        # Create a NumPy array
        Graphene = np.zeros([len(df), 7], dtype=object)
        
        for i in range(7):
            for j in range(len(df)):
                # Use the correct DataFrame column names or indices
                Graphene[j][i] = df.iloc[j, i]  # Assuming i is the correct column index
        
        # Convert the first five columns to integers
        Graphene[:, :5] = Graphene[:, :5].astype(int)
        
        # Save the NumPy array to a file
        np.savetxt("OLD_dihedral.itp", Graphene, delimiter=" ", fmt="%8d %6d %4d %7s %7s %7s %7s")



        #print(df)

        


        # Convert the 'ai' and 'aj' columns to string type if they're not already
        df['ai'] = df['ai'].astype(str)
        df['aj'] = df['aj'].astype(str)
        df['ak'] = df['ak'].astype(str)
        df['ah'] = df['ah'].astype(str)
        
        # Update 'ai' and 'aj', 'ak' columns using the exchange_dict
        df['ai'] = df['ai'].map(self.exchange_dict).fillna(df['ai'])
        df['aj'] = df['aj'].map(self.exchange_dict).fillna(df['aj'])
        df['ak'] = df['ak'].map(self.exchange_dict).fillna(df['ak'])
        df['ah'] = df['ah'].map(self.exchange_dict).fillna(df['ah'])

        #print(df)
        
        # Create a NumPy array
        Graphene = np.zeros([len(df), 7], dtype=object)
        
        for i in range(7):
            for j in range(len(df)):
                # Use the correct DataFrame column names or indices
                Graphene[j][i] = df.iloc[j, i]  # Assuming i is the correct column index
        
        # Convert the first five columns to integers
        Graphene[:, :5] = Graphene[:, :5].astype(int)
        
        # Save the NumPy array to a file
        np.savetxt("New_dihedral.itp", Graphene, delimiter=" ", fmt="%8d %6d %4d %7s %7s %7s %7s")
        
        return df




if __name__ == '__main__':
    itp = Itp("Protein_A.itp")
    #print("Atom SECTION:\n")
    #print(itp.atoms_extra)
    #print("BOND SECTION:\n")
    #print(itp.bonds)
    #print("ANGLE SECTION:\n")
    #print(itp.angles)
    #print("DIHEDRAL SECTION:\n")
    #print(itp.dihedrals)


### Box info 

file1 = open("header_Atom.txt","w")
file1.write("\n[ moleculetype ] \n") 
file1.write(" Protein_A         1 \n") 
file1.write("    ") 
file1.write("    ")
file1.write("    ")
file1.write("\n [ atoms ] \n") 
file1.close() 


file1 = open("header_Bond.txt","w")
file1.write("\n[ bonds ] \n") 
file1.write(" ; Backbone bonds \n") 
file1.close() 


file1 = open("header_constraints.txt","w")
file1.write("\n[ constraints ] \n") 
file1.close() 


file1 = open("header_Angle.txt","w")
file1.write("\n[ angles ] \n") 
file1.write(" ; ; Backbone angles \n") 
file1.close() 



file1 = open("header_Dihedral.txt","w")
file1.write("\n[ dihedrals ] \n") 
file1.write(" ; ; Backbone dihedrals \n") 
file1.close() 


filenames = ["header_Atom.txt", "New_Atom.itp"  ,"header_Bond.txt", "New_bond.itp" , "header_constraints.txt" , "New_constraint.itp" ,  "header_Angle.txt", "New_angle.itp" ,"header_Dihedral.txt", "New_dihedral.itp"   ]

with open("Topology.itp", "w") as outfile:
    for filename in filenames:
        with open(filename) as infile:
            contents = infile.read()
            outfile.write(contents)
            
            

print("End of the code . good luck")

