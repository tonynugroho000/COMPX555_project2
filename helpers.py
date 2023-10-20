# ------------------- HELPER FUNCTIONS ---------------------
def get_active_aa_count(protein_sequence):
    # Define a dictionary of active site amino acids
    active_site_aa = set(['C', 'H', 'D', 'K', 'S', 'T', 'Y'])
    
    # Initialize a count to keep track of active site amino acids
    count = 0
    
    # Iterate through each character in the protein sequence
    for aa in protein_sequence:
        if aa in active_site_aa:
            count += 1
    
    return count

def get_binding_aa_count(protein_sequence):
    # Define a set of binding site amino acids
    binding_site_aa = set(['E', 'G', 'D'])
    
    # Initialize a count to keep track of binding site amino acids
    count = 0
    
    # Iterate through each character in the protein sequence
    for aa in protein_sequence:
        if aa in binding_site_aa:
            count += 1
    
    return count

# This is a simplified function to calculate the net PK of an amino acid chain
# Assumes it is operating in a neutral environment ph = 7 (as normal bodily PH is around 7.35 to 7.45)
def calculate_protein_pK(protein_sequence, pH):
    # Define pK values for relevant amino acid side chains
    pK_values = {
        'D': 3.9,  # Aspartic acid
        'E': 4.2,  # Glutamic acid
        'K': 10.5,  # Lysine
        'R': 12.5,  # Arginine
        'H': 6.0  # Histidine
    }
    
    # Initialize variables to keep track of protonated and deprotonated forms
    protonated = 0
    deprotonated = 0

    # Iterate through the amino acid sequence
    for amino_acid in protein_sequence:
        if amino_acid in pK_values:
            # Determine the charge of the amino acid side chain at the given pH
            if pH < pK_values[amino_acid]:
                protonated += 1
            else:
                deprotonated += 1
    
    # Calculate the net charge
    net_charge = protonated - deprotonated
    
    return net_charge

def longest_repeat_length(protein_sequence):
    max_repeat_length = 1 # holds the highest current number
    current = 1 # stores the current string being observed

    for i in range(1, len(protein_sequence)):
        if (protein_sequence[i] == protein_sequence[i-1]):
            current += 1
            if (current > max_repeat_length):
                max_repeat_length = current
        else:
            current = 1

    return max_repeat_length

def count_hydrogen_donors(protein_sequence):
    # Define amino acids with hydrogen donor atoms in their side chains
    donor_amino_acids = set(['R', 'K', 'T'])
    
    # Initialize the count
    count = 0
    
    # Iterate through the amino acid sequence and count donors
    for amino_acid in protein_sequence:
        if amino_acid in donor_amino_acids:
            count += 1
    
    return count

def count_hydrogen_acceptors(protein_sequence):
    # Define amino acids with hydrogen acceptor atoms in their side chains
    acceptor_amino_acids = set(['D', 'E'])
    
    # Initialize the count
    count = 0
    
    # Iterate through the amino acid sequence and count acceptors
    for amino_acid in protein_sequence:
        if amino_acid in acceptor_amino_acids:
            count += 1
    
    return count



# -------------------- MAIN FUNCTION -----------------------
def main():
    # Get active amino acid count for demo amino acid
    demo_acid = "WCDYRLRWDPRDYEGLWVLR"
    ph = 7;
    print(f"For the amino acid {demo_acid} -------------------")
    print("Number of AA involved in active activity: ", get_active_aa_count(demo_acid))
    print("Number of AA involved in binding", get_binding_aa_count(demo_acid))
    print("The simplified PK of an amino acid chain", calculate_protein_pK(demo_acid, ph))
    print("Longest repeat length: ", longest_repeat_length(demo_acid))
    print(f"Number of hydrogen donors: {count_hydrogen_donors(demo_acid)}")
    print(f"Number of hydrogen acceptors: {count_hydrogen_acceptors(demo_acid)} \n")


    demo_acid_2 = "KGPELGLSQFCGSLKQAAPA -----------------"
    print(f"For the amino acid {demo_acid_2}")
    print("Number of AA involved in active activity: ", get_active_aa_count(demo_acid_2))
    print("Number of AA involved in binding", get_binding_aa_count(demo_acid_2))
    print("The simplified PK of an amino acid chain", calculate_protein_pK(demo_acid_2, ph))
    print("Longest repeat length: ", longest_repeat_length(demo_acid_2))
    print(f"Number of hydrogen donors: {count_hydrogen_donors(demo_acid_2)}")
    print(f"Number of hydrogen acceptors: {count_hydrogen_acceptors(demo_acid_2)}")
    return

# ---------------------- EXECUTION -------------------------

if __name__ == "__main__":
    main()
