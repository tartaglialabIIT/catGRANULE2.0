import os
from catgranuleFunctions import *


c_o_c = ['Normalizedfrequencyofbeta-sheetinall-betaclass-Palauetal--Int-J-PeptideProteinRes-1981-19-394-401',
 'Aggregationmed-Tartaglia-JMolBiol2008-380-2--425-36',
 'Normalizedfrequencyofcoil-Nagano-J-Mol-Biol-1973-75-401-420',
 'Aggregationhigh-Ferdandez-Escamilla-Nat-Biotechnol-2004-22-10--1302-6',
 'Membrane-buriedpreferenceparameters-Argosetal--Eur-J-Biochem-1982-128-565-575',
 'Hydrophobicity-Rao_Argos-Biochim-Biophys-Acta1986-869-197-214',
 'Averagerelativeprobabilityofbeta-sheet-Kanehisa-Tsong-Biopolymers1980-19-1617-1628',
 'Normalizedfrequencyofbeta-sheet-Chou-Fasman-Adv-Enzymol-1978-47-45-148',
 'Normalizedfrequencyofmiddlehelix-Crawfordetal--Proc-Natl-Acad-Sci-USA1973-70-538-542',
 'Compositionofaminoacidsinmembraneproteins-percent--Cedanoetal--J-Mol-Biol-1997-266-594-600',
 'Normalizedfrequencyofbeta-sheet-Crawfordetal-Proc-Natl-Acad-Sci-USA1973-70-538-542',
 'Averagevolumeofburiedresidue-Chothia-Nature1975-254-304-308',
 'Meanareaburiedontransfer-Roseetal--Science1985-229-834-838',
 'Hydrophobicity-Sweetetal--J-Mol-Biol-1983-171-479-488',
 'NucleicAcidBinding-interface_center-Terribilinietal--RNA2006-12-1450-1462',
 'Relativefrequencyinalpha-helix-Prabhakaran-Biochem-J-1990-269-691-696',
 'Normalizedfrequencyofbeta-sheet-withweights-Levitt-Biochemistry1978-17-4277-4285',
 'charge',
 'Aggregation-Tartaglia-J-Mol-Biol-2010-402-919',
 'Transmembraneregionsofnon-mt-proteins-Nakashimaetal--Proteins1990-8-173-178',
 'Averagedturnpropensitiesinatransmembranehelix-Monneetal--J-Mol-Biol-1999-293-807-814',
 'fg',
 'Informationmeasureforpleated-sheet-Robson-Suzuki-J-Mol-Biol-1976-107-327-356',
 'NucleicAcidBinding-classicalRBD-Castelloetal--Cell2011-149-1393-1406',
 'Normalizedfrequencyofbeta-sheet-unweighted-Levitt-Biochemistry1978-17-4277-4285',
 'Normalizedfrequencyofcoil-Tanaka-Scheraga-Macromolecules1977-10-9-20',
 'Knowledge-basedmembrane-propensityscalefrom3D_HelixinMPtopodatabases-Punta-Maritan-Proteins2003-50-114-121',
 'Ratioofburiedandaccessiblemolarfractions-Janin-Nature1979-277-491-492',
 'Hydrophobicity-Black-Anal-Biochem-1991-193-72-82',
 'Beta-sheet-Deleage-Roux-ProteinEngineering1987-1-289-294',
 'Energytransferfromouttoin-95_buried--Radzicka-Wolfenden-Biochemistry1988-27-1664-1670',
 'TOP-IDB-DunkerAK-ProteinPeptLett-2008-15-9--956',
 'Aggregation-Pawar-J-Mol-Biol-2005-350-379',
 'Normalizedcompositionofmembraneproteins-Nakashimaetal--Proteins1990-8-173-178',
 'Hydrophobicity-Kyte_Doolittle-J-Mol-Biol-1982-157-105-132',
 'Hydrophobicity-Eisenbergetal-J-Mol-Biol-1984-179-125-142',
 'Hydrophobicity-Fauchereetal--Eur-J-Med-Chem-1983-18-369-375',
 'Aggregation-Conchillo-Sole-BMCBioinformatics2007-8-65',
 'Aggregation-Tartaglia-ProteinScience2005-14-2735-2740',
 'UnfoldOverFold-DunkerAK-ProteinPeptLett-2008-15-9--956',
 'Aggregationlow-Ferdandez-Escamilla-Nat-Biotechnol-2004-22-10--1302-6',
 'Hydrophobicity-Janin-Nature1979-277-491-492',
 'Turnpropensityscalefortransmembranehelices-Monneetal--J-Mol-Biol-1999-288-141-145',
 'Normalizedfrequencyofalpha-helix-Burgessetal-Isr-J-Chem-1974-12-239-286',
 'B-Value-DunkerAK-ProteinPeptLett-2008-15-9--956',
 'TheChou-Fasmanparameterofthecoilconformation-Charton-Charton-J-Theor-Biol-1983-111-447-450',
 'Normalizedfrequencyofalpha-helixinall-alphaclass-Palauetal--Int-J-PeptideProteinRes-1981-19-394-401',
 'Percentageofburiedresidues-Janinetal--J-Mol-Biol-1978-125-357-386',
 'NucleicAcidBinding-_interface_close-1-Terribilinietal--RNA2006-12-1450-1462',
 'Informationmeasureforalpha-helix-Robson-Suzuki-J-Mol-Biol-1976-107-327-356',
 'NucleicAcidBinding-interface-NucleicAcidsRes-2011-39-D277',
 'Relativefrequencyinbeta-sheet-Prabhakaran-Biochem-J-1990-269-691-696',
 'Aggregationhigh-Tartaglia-JMolBiol2008-380-2--425-36',
 'Hydrophobicity-Abraham_Leo-Proteins-Structure-FunctionandGenetics1987-2-130-152',
 'rg',
 'Normalizedrelativefrequencyofcoil-Isogaietal--Biopolymers1980-19-1183-1210',
 'NucleicAcidBinding-interface_close+1-Terribilinietal--RNA2006-12-1450-1462',
 'Averageflexibilityindices-Bhaskaran-Ponnuswamy-Int-J-PeptideProteinRes-1988-32-241-255',
 'Coil-Deleage-Roux-ProteinEngineering1987-1-289-294',
 'Aggregationlow-Tartaglia-JMolBiol2008-380-2--425-36',
 'Proportionofresidues95_buried-Chothia-J-Mol-Biol-1976-105-1-14',
 'Normalizedfrequencyofalpha-helix-unweightedLevitt-Biochemistry1978-17-4277-4285',
 'NucleicAcidBinding-nonclassicalRBD-Castelloetal--Cell2011-149-1393-1406',
 'Normalizedfrequencyofalpha-helix-Chou-Fasman-Adv-Enzym-1978-47-45-148',
 'Aggregationmed-Ferdandez-Escamilla-Nat-Biotechnol-2004-22-10--1302-6',
 'Averagerelativeprobabilityofinnerbeta-sheet-Kanehisa-Tsong-Biopolymers1980-19-1617-1628',
 'Propensitytobeburiedinside-Wertz-Scheraga-Macromolecules1978-11-9-15',
 'DisProt-DunkerAK-ProteinPeptLett-2008-15-9--956',
 'Hydrophobicity-Roseman-J-Mol-Biol-1988-200-513-522',
 'Alpha-helix-Deleage-Roux-ProteinEngineering1987-1-289-294',
 'Meanvolumesofresiduesburiedinproteininteriors-Harpazetal--Structure1994-2-641-649',
 'Averagerelativeprobabilityofhelix-Kanehisa-Tsong-Biopolymers1980-19-1617-1628',
 'NucleicAcidBinding-mRNAinteractome-Castelloetal--Cell2011-149-1393-1406',
 'NucleicAcidBinding-HOH-NucleicAcidsRes-2011-39-D277',
 'NucleicAcidBinding-HB-NucleicAcidsRes-2011-39-D277',
 'Normalizedfrequencyofalpha-helix-withweights-Levitt-Biochemistry1978-17-4277-4285',
 'NucleicAcidBinding-unknownRBD-Castelloetal--Cell2011-149-1393-1406',
 'Hydrophobicity-Bull_Breese-Arch-Biochem-Biophys-1974-161-665-670',
 'Transmembraneregionsofmt-proteins-Nakashimaetal--Proteins1990-8-173-178',
 'AAcompositionofmembraneproteins-Nakashimaetal--Proteins1990-8-173-178',
 'Knowledge-basedmembrane-propensityscalefrom1D_HelixinMPtopodatabases-Punta-Maritan-Proteins2003-50-114-121',
 'Proportionofresidues100_buried-Chothia-J-Mol-Biol-1976-105-1-14',
 'n_contacts',
 'RG_protein',
 'rmsd_prot',
 'Length',
 'fullCharge',
 'extCharge',
 'Percentage_Coil_FullSeq',
 'Percentage_AlphaHelix_FullSeq',
 'Percentage_BetaBridge_FullSeq',
 'Percentage_Strand_FullSeq',
 'Percentage_3_10Helix_FullSeq',
 'Percentage_PiHelix_FullSeq',
 'Percentage_Turn_FullSeq',
 'Percentage_Bend_FullSeq',
 'Percentage_Coil_ExtSeq',
 'Percentage_AlphaHelix_ExtSeq',
 'Percentage_BetaBridge_ExtSeq',
 'Percentage_Strand_ExtSeq',
 'Percentage_3_10Helix_ExtSeq',
 'Percentage_PiHelix_ExtSeq',
 'Percentage_Turn_ExtSeq',
 'Percentage_Bend_ExtSeq',
 'asa_mean',
 'asa_std',
 'RBD_Full_min_0',
 'RBD_Full_min_1',
 'RBD_Full_min_2',
 'RBD_Full_indexmin_0',
 'RBD_Full_indexmin_1',
 'RBD_Full_indexmin_2',
 'RBD_int_min_0',
 'RBD_int_min_1',
 'RBD_int_min_2',
 'RBD_int_indexmin_0',
 'RBD_int_indexmin_1',
 'RBD_int_indexmin_2',
 'RBD_ext_min_0',
 'RBD_ext_min_1',
 'RBD_ext_min_2',
 'RBD_ext_indexmin_0',
 'RBD_ext_indexmin_1',
 'RBD_ext_indexmin_2',
 'average_plddt',
 'stddev_plddt',
 'n_contacts_norm',
 'RG_protein_norm']
correct_order_columns = np.array(c_o_c, dtype='str')


def compute_score_and_profile_from_aa_text(aa_text):
    '''Computes the Profiles of the sequences and the predictions of LLPS
    Input: an aminoacid string string
    Output: a Profiles array, a Predictions array'''
    code_dir = os.getcwd()
    scales_dir = './src/ChemicalPhysicalScales_Py_dictionary'
    classifiers_dir = './src/TRAINED_MODELS/'
    # get the sequences and the according ids
    sequences, ids = get_sequences_from_aa_text(aa_text)
    # get the dataframe with the physical/chemical properties
    pc_df = get_physical_chemical_properties(sequences, ids, scales_dir)
    # reorder the columns with the correct order
    pc_df = pc_df[correct_order_columns[:82]]
    # get the matrixes needed in order to compute the profiles
    profiles = get_physical_chemical_profiles(sequences, scales_dir, classifiers_dir, correct_order_columns)
    # get the llps predictions
    predictions = predict(pc_df, classifiers_dir, only_pc=True)
    prof = [(smooth(profiles[i],21)/np.max((smooth(profiles[i],21))))*predictions['RandomForest'][i] for i in range(len(profiles))]

    # prof = [(profiles[i]/np.max(profiles[i]))*predictions[i] for i in range(len(profiles))]
    # create dataframes
    #profiles_df = pd.DataFrame(data=profiles, index=ids)
    profiles_dict = dict(zip(ids, prof))
    predictions_df = pd.DataFrame(data=predictions, index=ids)
    # return the according arrays
    return profiles_dict, predictions_df


def compute_score_and_profile_from_text(fasta_text):
    '''Computes the Profiles of the sequences and the predictions of LLPS
    Input: a fasta string
    Output: a Profiles array, a Predictions array'''
    code_dir = os.getcwd()
    scales_dir = './src/ChemicalPhysicalScales_Py_dictionary'
    classifiers_dir = './src/TRAINED_MODELS/'
    # get the sequences and the according ids
    sequences, ids = get_sequences_from_fasta(fasta_text)
    # get the dataframe with the physical/chemical properties
    pc_df = get_physical_chemical_properties(sequences, ids, scales_dir)
    # reorder the columns with the correct order
    pc_df = pc_df[correct_order_columns[:82]]
    # get the matrixes needed in order to compute the profiles
    profiles = get_physical_chemical_profiles(sequences, scales_dir, classifiers_dir, correct_order_columns)
    # get the llps predictions
    predictions = predict(pc_df, classifiers_dir, only_pc=True)
    # create dataframes
    prof = [(smooth(profiles[i],21)/np.max((smooth(profiles[i],21))))*predictions['RandomForest'][i] for i in range(len(profiles))]
    profiles_dict = dict(zip(ids, prof))
    predictions_df = pd.DataFrame(data=predictions, index=ids)
    # return the according arrays
    return profiles_dict, predictions_df


def compute_score_and_profile_from_fasta_file(fasta_file):
    '''Computes the Profiles of the sequences and the predictions of LLPS
    Input: a fasta file
    Output: a Profiles array, a Predictions array'''
    code_dir = os.getcwd()
    scales_dir = './src/ChemicalPhysicalScales_Py_dictionary'
    classifiers_dir = './src/TRAINED_MODELS/'
    # get the sequences and the according ids
    sequences, ids = get_sequences_from_file(fasta_file)
    # get the dataframe with the physical/chemical properties
    pc_df = get_physical_chemical_properties(sequences, ids, scales_dir)
    # reorder the columns with the correct order
    pc_df = pc_df[correct_order_columns[:82]]
    # get the matrixes needed in order to compute the profiles
    profiles = get_physical_chemical_profiles(sequences, scales_dir, classifiers_dir, correct_order_columns)
    # get the llps predictions
    predictions = predict(pc_df, classifiers_dir, only_pc=True)
    # create dataframes
    prof = [(smooth(profiles[i],21)/np.max((smooth(profiles[i],21))))*predictions['RandomForest'][i] for i in range(len(profiles))]
    profiles_dict = dict(zip(ids, prof))
    predictions_dict = dict(zip(ids, predictions))
    predictions_df = pd.DataFrame(data=predictions, index=ids)
    # return the according arrays
    # return profiles_dict, predictions_df
    return prof, predictions['RandomForest']

def compute_score_and_profile_from_fasta_file_for_mutscan(fasta_file):
    '''Computes the Profiles of the sequences and the predictions of LLPS
    Input: a fasta file
    Output: a Profiles array, a Predictions array'''
    code_dir = os.getcwd()
    scales_dir = './src/ChemicalPhysicalScales_Py_dictionary'
    classifiers_dir = './src/TRAINED_MODELS/'
    # get the sequences and the according ids
    sequences, ids = get_sequences_from_file(fasta_file)
    # get the dataframe with the physical/chemical properties
    pc_df = get_physical_chemical_properties(sequences, ids, scales_dir)
    # reorder the columns with the correct order
    pc_df = pc_df[correct_order_columns[:82]]
    # get the matrixes needed in order to compute the profiles
    profiles = get_physical_chemical_profiles(sequences, scales_dir, classifiers_dir, correct_order_columns)
    # get the llps predictions
    predictions = predict(pc_df, classifiers_dir, only_pc=True)
    # create dataframes
    #profiles_dict = dict(zip(ids, prof))
    #predictions_dict = dict(zip(ids, predictions))
    #predictions_df = pd.DataFrame(data=predictions, index=ids)
    # return the according arrays
    # return profiles_dict, predictions_df
    return profiles, predictions['RandomForest']


def compute_score_and_profile_from_pdb(pdb_files_string):
    '''Computes the Profiles of the sequences and the predictions of LLPS
    Input: a directory with pdb files
    Output: a Profiles array, a Predictions array'''
    code_dir = os.getcwd()
    scales_dir = './src/ChemicalPhysicalScales_Py_dictionary'
    classifiers_dir = './src/TRAINED_MODELS/'
    af_code_dir = './src/AlphaFold/'
    # convert the pdb string with the files in a list of paths
    pdb_files = get_list_from_pdb_string(pdb_files_string)
    # get the sequences and the according ids
    sequences, ids = get_sequences_from_pdbs(pdb_files)
    # get the dataframe with the physical/chemical properties
    pc_df = get_physical_chemical_properties(sequences, ids, scales_dir)
    # get the correct index | UNCOMMENT THIS IF YOU HAVE AF MODELS
    #pc_df.index = pc_df.index.map(return_af_protein_name)
    # get the dataframe with the alphafold properties
    af_df = get_alphafold_features_from_file(af_code_dir, pdb_files)
    # add plddt columns
    af_complete_df = add_plddt_to_af_features_from_file(pdb_files, af_df)
    # postprocess alphafold features
    postprocessed_af_df = ProcessAF_Features_from_file(af_complete_df)
    # set the uniprot _ID as the index
    postprocessed_af_df.set_index('Uniprot_ID', inplace=True)
    # get the matrixes needed in order to compute the profiles
    profiles = get_physical_chemical_profiles(sequences, scales_dir, classifiers_dir, correct_order_columns)
    # create the combined dataframe of the physical/chemical properties
    cmb_df = pd.concat([pc_df, postprocessed_af_df], axis=1)
    # reorder the columns in order to get the correct columns
    cmb_df = cmb_df[correct_order_columns]
    # get the llps predictions
    predictions = predict(cmb_df, classifiers_dir)
    # create dataframes
    prof = [(smooth(profiles[i],21)/np.max((smooth(profiles[i],21))))*predictions['MLP'][i] for i in range(len(profiles))]

    profiles_dict = dict(zip(ids, prof))
    predictions_df = pd.DataFrame(data=predictions, index=ids)

    profiles_df=pd.DataFrame([profiles_dict]).transpose()
    scoPd = predictions_df.reset_index()
    proPd = profiles_df.reset_index()
    scoPd.columns=['FileName', 'LLPS_Score']
    proPd.columns=['FileName', 'LLPS_Profile']
    # return the according arrays
    finalDf = pd.merge(scoPd, proPd, on='FileName')
    return finalDf


def compute_score_profile_ListPdb(pdbList):
    sc=[]
    pr=[]
    nm=[]
    for pdName in pdbList:

        df =compute_score_and_profile_from_pdb(pdName)
        nm.append(df.FileName[0])
        sc.append(df.LLPS_Score[0])
        pr.append(df.LLPS_Profile[0])

    finalDf = pd.DataFrame({
        'FileName': nm,
        'LLPS_Score': sc,
        'LLPS_Profile': pr,
    })

    return finalDf



def compute_score_profile_fatsa_DF(fasta_file):

    prof, scores  = compute_score_and_profile_from_fasta_file(fasta_file)
    names = [i.id for i in list(SeqIO.parse(fasta_file, "fasta"))]
    
    finalDf = pd.DataFrame({
        'Name': names,
        'LLPS_Score': scores,
        'LLPS_Profile': prof,
    })

    return finalDf

def compute_score_profile_fatsa_DF_for_mutscan(fasta_file):

    prof, scores  = compute_score_and_profile_from_fasta_file_for_mutscan(fasta_file)
    names = [i.id for i in list(SeqIO.parse(fasta_file, "fasta"))]
    
    finalDf = pd.DataFrame({
        'Name': names,
        'LLPS_Score': scores,
        'LLPS_Profile': prof,
    })

    return finalDf


def catGranule2_str(aa_text):
    '''Computes the Profiles of the sequences and the predictions of LLPS
    Input: an aminoacid string string
    Output: a Profiles array, a Predictions array'''
    code_dir = os.getcwd()
    scales_dir = './src/ChemicalPhysicalScales_Py_dictionary'
    classifiers_dir = './src/TRAINED_MODELS/'
    # get the sequences and the according ids
    sequences, ids = get_sequences_from_aa_text(aa_text)
    # get the dataframe with the physical/chemical properties
    pc_df = get_physical_chemical_properties(sequences, ids, scales_dir)
    # reorder the columns with the correct order
    pc_df = pc_df[correct_order_columns[:82]]
    # get the matrixes needed in order to compute the profiles
    profiles = get_physical_chemical_profiles(sequences, scales_dir, classifiers_dir, correct_order_columns)[0]
    # get the llps predictions
    predictions = predict(pc_df, classifiers_dir, only_pc=True)['RandomForest'][0]
    prof = [(smooth(profiles[i],21)/np.max((smooth(profiles[i],21))))*predictions['RandomForest'][i] for i in range(len(profiles))]

    return predictions, prof
