"""Program converts helices occuring in protein to ones containing
fragment with beta-amino acids, according to pattern. Program can be imported
as module.

File contains following functions:

    * create_pattern_set - creates all possible variants of pattern sequence
    * which_pattern_set - determines number of pattern in the set
    * generate_new_sequence - generates sequences with non-standard amino acids
    * write_torsion_angles - saving values of torsion angles of original protein
    * set_standard_angle - sets values of torsion angles for standard residues
    * set_nonstandard_angle - sets values of torsion angles for non-standard
                              residues
    * save_to_file - saves (.pdb) file with mutated protein
    * swap_helix - swaps particular helix with fragment containing non-standard
                   amino acids
    * swap_all - swaps all of the helices with fragment containing non-standard
                   amino acids
    * main - the main function of the script
"""

from copy import deepcopy
from pyrosetta import *
from pyrosetta.toolbox import cleanATOM
init()

def create_pattern_set():
    """Generates variuos variants of pattern.

    Returns
    -------
    pattern_set : list
        list of different variant of pattern sequence
    torsion_angles_set : list
        list of torsion angles corrensponding to particular variants
    """

    pattern_set = []
    torsion_angles_set = []

    pattern = "AAXAAAX"
    pattern_list = list(pattern)
    torsion_angles = [[-62.29, -47.18, -176.91, 0.0],
            [-64.66, -48.01, 180.0, 0.0],
            [-104.04, -110.08, -179.96, 83.13],
            [-59.21, -44.27, 178.53, 0.0],
            [-65.47, -44.49, -179.40, 0.0],
            [-63.85, -42.24, 179.99, 0.0],
            [-118.65,-103.67, -179.96, 83.46]]

    for aa in range(len(pattern_list)):
        #creating a list of all possible variants of pattern
        pattern_list.append(pattern_list.pop(0))
        pattern_set.append(''.join(pattern_list))

        #creating a list of tosion angles corresponding to particular variants
        torsion_angles.append(torsion_angles.pop(0))
        torsion_angles_copy = deepcopy(torsion_angles)
        torsion_angles_set.append(torsion_angles_copy)

    return pattern_set, torsion_angles_set

def which_pattern_set(sequence, pattern_set):
    """Determines number of set for pattern and torsion angles.

    Parameters
    ----------
    sequence : str
        sequence of fragment with beta amino acids passed by user
    pattern_set : list
        list of different variant of pattern sequence

    Returns
    -------
    set_number : int
        number of set
    """
    for i in range(len(pattern_set)):
        if(sequence[0:7] == pattern_set[i]):
            set_number = i
            return set_number

def generate_new_sequence(sequence):
    """Generates sequence with beta amino acids which is readable for pyrosetta
    software. X -> A[transACPC]

    Parameters
    ----------
    sequence : str
        sequence of fragment with beta amino acids passed by user

    Returns
    -------
    pyrosetta_sequence : list
        sequence on beta amino acids readable for pyrosetta
    """
    pyrosetta_sequence = []
    for i in sequence:
        if(i == 'X'):
            i = 'A[transACPC]'
            pyrosetta_sequence.append(i)
        else:
            pyrosetta_sequence.append(i)

    return pyrosetta_sequence

def write_torsion_angles(protein):
    """Collects values of torsion angles for original protein.

    Parameters
    ----------
    protein : Pose (PyRosetta class)

    Returns
    -------
    phi, psi, omega : list
        lists of phi, psi and omega torsion angles
    """
    phi = []
    psi = []
    omega = []

    for i in range(1, protein.total_residue() + 1):
        phi.append(protein.phi(i))
        psi.append(protein.psi(i))
        omega.append(protein.omega(i))

    return phi, psi, omega

def set_standard_angle(pose_mutant, i, phi, psi, omega):
    """Sets values of torsion angles.

    Parameters
    ----------
    pose_mutant : Pose (PyRosetta class)
        object of mutated protein
    i : int
        control variable
    phi, psi, omega : list
        lists of phi, psi and omega torsion angles

    Returns
    -------
    pose_mutant : Pose (PyRosetta class)
        object of mutated protein
    """
    pose_mutant.set_phi(i, phi[i-1])
    pose_mutant.set_psi(i, psi[i-1])
    pose_mutant.set_omega(i, omega[i-1])
    return pose_mutant

def set_nonstandard_angle(pose_mutant, i, counter, angles, set_number):
    """Sets values of torsion angles for nonstandard residues.

    Parameters
    ----------
    pose_mutant : Pose (PyRosetta class)
        object of mutated protein
    i : int
        control variable
    counter : int
        control variable
    angles : list
        list of torsion angles for nonstandard residues
    set_number : int
        number of set

    Returns
    -------
    pose_mutant : Pose (PyRosetta class)
        object of mutated protein
    """
    pose_mutant.set_phi(i, angles[set_number][(counter)%7][0])
    pose_mutant.set_psi(i, angles[set_number][(counter)%7][1])
    pose_mutant.set_omega(i, angles[set_number][(counter)%7][2])
    return pose_mutant

def save_to_file(pose_mutant, helix):
    """Saves (.pdb) file with mutated protein.

    Parameters
    ----------
    pose_mutant : Pose (PyRosetta class)
        object of mutated protein
    helix : int
        number of helix
    """
    fileToSave = "mutant_" + str(helix) + ".pdb"
    #Writing data to (.pdb) file
    pose_mutant.dump_pdb(fileToSave)

def swap_helix(protein, borders, protein_seq, pyrosetta_seq, phi, psi, omega,
               angles, set_number, fragment):
    """Swaps particular helix.

    Parameters
    ----------
    protein : Pose (PyRosetta class)
    borders : list
        list of end points for 3 helices
    protein_seq : str
        sequence of original protein
    pyrosetta_seq : list
        sequence of fragment to be inserted
    phi, psi, omega : list
        lists of phi, psi and omega torsion angles
    angles : list
        list of torsion angles for nonstandard residues
    set_number : int
        number of set
    """
    #Creating list for sequence of mutated protein
    mutant = []
    counter = 0

    #Determination of sequence
    for i in range(1, protein.total_residue()+1):
        #Fragment with standard amino acids
        if (i < borders[fragment - 1][0]) or (i > borders[fragment - 1][1]):
            mutant.append(protein_seq[i-1])
        #Fragment with nonstandard amino acids
        else:
            mutant.append(pyrosetta_seq[counter%7])
            counter += 1
    #Saving data to pose object
    pose_mutant = pose_from_sequence(''.join(mutant), 'fa_standard')

    counter = 0
    #Determination of torsion angles
    for i in range(1, protein.total_residue() + 1):
        #Fragment with standard amino acids
        if (i < borders[fragment - 1][0]) or (i > borders[fragment - 1][1]):
            pose_mutant = set_standard_angle(pose_mutant, i, phi, psi, omega)
        #Fragment with nonstandard amino acids
        else:
            pose_mutant = set_nonstandard_angle(pose_mutant, i, counter,
                                                angles, set_number)
            #Setting value for additional beta amino acid angle
            if (pose_mutant.residue(i).name() == 'transACPC'):
                pose_mutant.set_theta(i, angles[set_number][(counter)%7][3])
            counter += 1

        save_to_file(pose_mutant, fragment)

def swap_all(protein, borders, protein_seq, pyrosetta_seq, phi, psi, omega,
             angles, set_number, fragment):
    """Swaps all helices in protein.

    Parameters
    ----------
    protein : Pose (PyRosetta class)
    borders : list
        list of end points for 3 helices
    protein_seq : str
        sequence of original protein
    pyrosetta_seq : list
        sequence of fragment to be inserted
    phi, psi, omega : list
        lists of phi, psi and omega torsion angles
    angles : list
        list of torsion angles for nonstandard residues
    set_number : int
        number of set
    """
    #Creating list for sequence of mutated protein
    mutant = []
    counter = 0

    #Determination of sequence
    for i in range(1, protein.total_residue()+1):
        #Fragment with standard amino acids
        if (i == 1 or i == 2 or (i >= 11 and i <= 13) or i == 20 or i == 21
           or (i >= 33 and i <= 35)):
            mutant.append(protein_seq[i-1])
        #Fragment with nonstandard amino acids
        else:
            mutant.append(pyrosetta_seq[counter%7])
            counter += 1
            if ((i==borders[0][-1]) or (i==borders[1][-1])
                 or (i==borders[2][-1])):
                counter = 0
    pose_mutant = pose_from_sequence(''.join(mutant), 'fa_standard')

    #Determination of torsion angles
    for i in range(1, protein.total_residue() + 1):
        #Fragment with standard amino acids
        if (i == 1 or i == 2 or (i >= 11 and i <= 13) or i == 20 or i == 21
            or (i >= 33 and i <= 35)):
            pose_mutant = set_standard_angle(pose_mutant, i, phi, psi, omega)
        #Fragment with nonstandard amino acids
        else:
            pose_mutant = set_nonstandard_angle(pose_mutant, i, counter,
                                                angles, set_number)
            if (pose_mutant.residue(i).name() == 'transACPC'):
                pose_mutant.set_theta(i, angles[set_number][(counter)%7][3])
            counter += 1

            if ((i==borders[0][-1]) or (i==borders[1][-1])
                 or (i==borders[2][-1])):
                counter = 0

        save_to_file(pose_mutant, fragment)

def main():

    #End points of protein helices
    border1N = 3
    border1C = 10

    border2N = 14
    border2C = 19

    border3N = 22
    border3C = 32

    borders = [[border1N, border1C], [border2N, border2C], [border3N, border3C]]

    pattern_set, torsion_angles_set = create_pattern_set()

    sequence = input("Insert sequence:")

    set_number = which_pattern_set(sequence, pattern_set)

    pyrosetta_sequence = generate_new_sequence(sequence)

    #Cleaning original protein (.pdb) file
    cleanATOM('1wy3_ww.pdb')
    pose = pose_from_sequence(''.join(pyrosetta_sequence), 'fa_standard')
    protein = pose_from_pdb('1wy3_ww.clean.pdb')

    phi, psi, omega = write_torsion_angles(protein)

    protein_sequence = list(protein.sequence())
    #Change in amino acid no. 23, NORleucine -> Leucine
    protein_sequence[23] = 'L'

    helix = int(input('Which helix do you want to swap?'))

    if(helix < 4):
        swap_helix(protein, borders, protein_sequence, pyrosetta_sequence,
                   phi, psi, omega, torsion_angles_set, set_number, helix)
    else:
        swap_all(protein, borders, protein_sequence, pyrosetta_sequence,
                 phi, psi, omega, torsion_angles_set, set_number, helix)

if __name__ == '__main__':
    main()
