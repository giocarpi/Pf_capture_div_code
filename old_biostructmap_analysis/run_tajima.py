import biostructmap

#Change fname for each protein and make sure to update pdb name too on line 8
fname = 'pfs230_seg13'
print(fname)

# Initialise structure object -- update pdb name between proteins
structure = biostructmap.Structure('input_data/pfs230_segment13_AF_cleaned.pdb', 'test_pdb_name')

# Read in multiple sequence alignment data 
msa_data = biostructmap.SequenceAlignment('input_data/'+fname+'_trimmed.fasta')
data = {('A',): msa_data}

# Reference seq -- 3d7 reference sequence is first sequence
reference_seq = {'A': str(msa_data[0].seq)}

# methods calculated are tajimasd and nucleotide diversity -- note rsa range goes above 1 as theoretical maximum of 1 can be exceeded
results1 = structure.map(data=data, method='tajimasd', ref=reference_seq, radius=15, map_to_dna=True, rsa_range=(0.15,99))
results2 = structure.map(data=data, method='nucleotide_diversity', ref=reference_seq, radius=15, map_to_dna=True, rsa_range=(0.15,99))

# save tajima results
results1.write_data_to_pdb_b_factor(fileobj='result_data/'+fname+'_trim_tajima.pdb')
results1.write_residue_data_to_csv(fileobj='result_data/'+fname+'_trim_tajima.csv')

# save watterman results
results2.write_data_to_pdb_b_factor(fileobj='result_data/'+fname+'_trim_pi.pdb')
results2.write_residue_data_to_csv(fileobj='result_data/'+fname+'_trim_pi.csv')