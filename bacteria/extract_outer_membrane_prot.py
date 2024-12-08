import os
from Bio import SeqIO
import subprocess

def parse_psortb_results(psortb_output_file):
    """Parse psortb results to extract SeqIDs with 'OuterMembrane' as Final Prediction."""
    outer_membrane_ids = []
    with open(psortb_output_file, 'r') as file:
        seq_id = None
        is_outer_membrane = False
        for line in file:
            if line.startswith("SeqID:"):
                if seq_id and is_outer_membrane:
                    outer_membrane_ids.append(seq_id)
                seq_id = line.split(" ", 1)[1].strip()
                is_outer_membrane = False  # Reset for the new sequence
            elif line.startswith("  Final Prediction:"):
                if "OuterMembrane" in line:
                    is_outer_membrane = True
        # Append the last seq_id if applicable
        if seq_id and is_outer_membrane:
            outer_membrane_ids.append(seq_id)
    return outer_membrane_ids


def filter_fasta_by_ids(fasta_file, output_file, ids_to_keep):
    """Filter a FASTA file to include only sequences with IDs in ids_to_keep."""
    with open(output_file, 'w') as out_f:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in ids_to_keep:
                SeqIO.write(record, out_f, "fasta")


def main(psortb_results_folder, fasta_folder, output_folder):
    """Main function to process multiple psortb result and proteome files."""
    os.makedirs('psortb_res', exist_ok=True)
    for fasta_file in os.listdir(fasta_folder):
        # 1) Generate Psortb file
        try:
            subprocess.run(f'./psortb_app -i {fasta_file} -r {psortb_results_folder} -n', shell = True, capture_output = True)
        except:
            print("Error: psortb command not executed")
            break
        
        # 2) Create a list of only membrane proteins ids
        for psortb_output_file in os.listdir(psortb_results_folder):
            outer_membrane_ids = parse_psortb_results(psortb_output_file)

        # 3) Create a new fasta file containing only membrane prot
        output_file = os.path.join(output_folder, f"{fasta_file}_only_OuterMembrane.fasta")
        filter_fasta_by_ids(fasta_file, output_file, outer_membrane_ids)
        print(f"Processed {fasta_file}: {len(outer_membrane_ids)} OuterMembrane proteins saved.")

if __name__ == "__main__":
    # Input folders
    PSORTB_RESULTS_FOLDER = "C:/Users/lorenzo/Desktop/results"
    FASTA_FOLDER = "C:/Users/lorenzo/Desktop/datasets/bacteria_proteomes"
    OUTPUT_FOLDER = "C:/Users/lorenzo/Desktop/path_to_fasta_files"
    
    main(PSORTB_RESULTS_FOLDER, FASTA_FOLDER, OUTPUT_FOLDER)

"""     for psortb_file in os.listdir(psortb_results_folder):
        if not psortb_file.endswith(".txt"):
            continue  # Skip non-result files
        bacterium_name = psortb_file.replace(".txt", "")
        psortb_path = os.path.join(psortb_results_folder, psortb_file)
        
        # Parse psortb results
        outer_membrane_ids = parse_psortb_results(psortb_path)
        
        # Find the corresponding FASTA file
        fasta_file = os.path.join(fasta_folder, f"{bacterium_name}.fasta")
        if not os.path.exists(fasta_file):
            print(f"FASTA file for {bacterium_name} not found. Skipping.")
            continue
        
        # Create output file
        output_file = os.path.join(output_folder, f"{bacterium_name}_only_OuterMembrane.fasta")
        filter_fasta_by_ids(fasta_file, output_file, outer_membrane_ids)
        print(f"Processed {bacterium_name}: {len(outer_membrane_ids)} OuterMembrane proteins saved.")
 """


