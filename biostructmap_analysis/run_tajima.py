import sys
import biostructmap

# check for five arguments (script name, fname, pdb_file, fasta_file, chain_id)
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print('\n\n')
        print("Usage: run_tajima.py <fout_name> <pdb_file> <fasta_file> <chain_id>")
        print('\n\n')
        sys.exit(1)
    else:
        _, fname, pdb_file, fasta_file, chain_id = sys.argv

#Change fname for each protein and make sure to update pdb name too on line 8
# fname = 'something'
print(fname)

# Initialise structure object -- update pdb name between proteins
pdb_file = 'pdbs/' + pdb_file
structure = biostructmap.Structure(pdb_file, 'test_pdb_name')

# Read in multiple sequence alignment data 
fasta_file = 'fastas/' + fasta_file
msa_data = biostructmap.SequenceAlignment(fasta_file)
data = {(chain_id,): msa_data}

# Reference seq -- 3d7 reference sequence is first sequence
reference_seq = {chain_id: str(msa_data[0].seq)}

# methods calculated are tajimasd and nucleotide diversity -- note rsa range goes above 1 as theoretical maximum of 1 can be exceeded
results1 = structure.map(data=data, method='tajimasd', ref=reference_seq, radius=15, map_to_dna=True, rsa_range=(0.15,99))
results2 = structure.map(data=data, method='nucleotide_diversity', ref=reference_seq, radius=15, map_to_dna=True, rsa_range=(0.15,99))

# save tajima results
results1.write_data_to_pdb_b_factor(fileobj='output/'+fname+'_trim_tajima.pdb')
results1.write_residue_data_to_csv(fileobj='output/'+fname+'_trim_tajima.csv')

# save watterman results
results2.write_data_to_pdb_b_factor(fileobj='output/'+fname+'_trim_pi.pdb')
results2.write_residue_data_to_csv(fileobj='output/'+fname+'_trim_pi.csv')
