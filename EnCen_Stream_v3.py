import streamlit as st
import numpy as np
import pandas as pd
import os
import os.path
import numpy as np 
import pandas as pd 
import numpy as np 
red = "\033[91m"
reset_color = "\033[0m"
import seaborn as sns
import matplotlib.pyplot as plt
from EnCen_Functions import diff_score, EC_extract, tsv_to_fasta, diamond_impl, genome_extractor_syn,  EC_corrector
import altair as alt
from streamlit.runtime.scriptrunner import add_script_run_ctx,get_script_run_ctx
from matplotlib.colors import Normalize
#_______________________________________________________________________________________
col1, col2, col3 = st.columns([1, 2, 1])
with col2: 
    st.image('/home/anna/Documents/EcoGenoRisk_Paper_Revisions/GUI_Output/logo2.png', width=800)
# st.markdown("<h1 style='text-align: center;'>EcoGenoRisk</h1>", unsafe_allow_html=True)
st.markdown("<h1 style='text-align: center; font-size: 23px; color: #0F4A0F; '>Metagenomic Synthetic Biology Threat Assessment</h1>", unsafe_allow_html=True)
# st.markdown("<h1 style='text-align: center; font-size: 25px; '>Developed by John Docter (john.docter@colorado.edu), University of Colorado - Boulder>", unsafe_allow_html=True)

# st.markdown('### Metagenomic Synthetic Biology Risk Assessment')
st.divider()
st.header('Instructions')
with st.expander("### Start Here"):
    st.write('This is a simple interface for the EcoGenoRisk pipeline.')
    st.write('User is required to have the amino acid (.faa) files of the synthetic and or comparison organism. Test files are available in the associated library for metagenome comparison')
    st.write('For larger or custom analyses, please refer to EcoGenoRisk Source Code :blue[https://github.com/UCBoulder/EcoGenoRisk]')
#Ask for home directory and where all the files should be saved______________________________________________________________________________________________________________
    st.divider()
    home_dir = st.text_input('Please Enter the filepath where you would like all outputs saved')
    # home_dir = '/home/anna/Documents/EcoGenoRisk_Paper_Revisions/GUI_Output'
    if not home_dir:
        st.stop()
    if '/' not in home_dir:
        st.write(':red[Is this a filepath?]')
    else:  
        st.write(':green[All outputs will be saved to ]  ')
        st.write(home_dir)



#Create Unitprot & EC_Library________________________________________________________________________________________
st.header('Enzyme and Protein Library Creation')
with st.expander('Data Creation'):
    if st.button('Press to create Expasy and Unitprot Libraries'):
        if os.path.isfile(home_dir +'/EC_library.csv'):
                st.success(":green[EC Library Detected]")
        else:
            with st.spinner('Creating EC_Library'):
                EC_extract()
            st.success('EC_Library Created')
            st.write(home_dir)

        if os.path.isfile(home_dir + '/uniprot.fasta'):
            st.success(':green[Protein Library Detected]')
        else: 
            with st.spinner('Creating Protein Reference Library'):
                tsv_to_fasta()
            st.success('Protein Reference Created')
    # else:
    #     st.stop()
    # else:
    #     st.stop()

#Upload Synbio .faa file____________________________________________________________________________
st.header('Synbio File Upload')
with st.expander('Upload'):
    uploaded_file_synbio = st.file_uploader(" ", key = 'IW_syn')
    if not uploaded_file_synbio: 
        st.stop()
    if uploaded_file_synbio:
        synbio_path = os.path.join(home_dir, uploaded_file_synbio.name)
        with open(synbio_path, 'wb') as f:
            f.write(uploaded_file_synbio.getvalue()) #uploaded file is in the home dir at this point 
    #Parsing to align with metagenomes______________________________________________________________________



    # #Synbio Diamond Processing and outputs synbio functional profile__________________________________________________________________________

    reference = (home_dir + '/uniprot.fasta')
    file_name = uploaded_file_synbio.name
    with st.spinner('Diamond Aligner Matching Synbio.faa Sequences to Unitprot Reference'):
        diamond_syn = diamond_impl(home_dir, file_name, reference) #returns the output folder, in this case home_dir




#Genome Extractor_______________________________________________________________________________________________
    output2 = genome_extractor_syn(diamond_syn, file_name, home_dir)
    st.success('Synbio Functional Profile Created')

    #Parsing to align with metagenomes______________________________________________________________________
    parsed_found = any('Parsed' in filename for filename in os.listdir(home_dir))
    if parsed_found: 
        st.success(':green[Parsed list found]')
    else:
        st.spinner('Parsing')
        EC_corrector(home_dir)


#Upload Pre-Made Metagenome Functional Profile and save to home directory_________________________________________________________
st.header('Metagenome Upload')
with st.expander('## Please choose from any of the pre-processed metagenomes. Multiples are ok'):
    
    file_found = True
    for all_files in os.listdir(home_dir):
        if 'Score' in all_files:
            st.success(':green[Scores Detected]')
            file_found = False
            break

    if file_found:
        uploaded_file_meta = st.file_uploader(" ", accept_multiple_files= True, key='IW')
        if not uploaded_file_meta:
            st.write('Waiting on File Upload')
            st.stop()

        if uploaded_file_meta:
            for every_file in uploaded_file_meta: 
                biome_path = os.path.join(home_dir, every_file.name)
                new_name = every_file.name + ' Difference Score'
                synbio_func_path = home_dir +'/' + 'BioPesticide Parsed.txt'
                with open(biome_path, 'wb') as f:
                    f.write(every_file.getvalue())
                with st.spinner('Working...'):
                    diff_score(synbio_func_path, biome_path)
                with open('Absolute_Difference_Comparison_Score.txt', 'r') as f:
                    file_content = f.read()
                os.rename('Absolute_Difference_Comparison_Score.txt', new_name )
    

    


#3) Then, header with Visualization for score. We'll just hide the fact that scores are being calculated in the upload. 
st.header('Scoring', divider='grey')
if st.button('Press to See Threat Assessment', key='IW_3'):
    top_ten_master_list = []
    for file in os.listdir(home_dir):
        if 'Score' in file:
            biome_name = "_".join(file.split("_")[:2])
            biome_name = biome_name.replace('_', " ")
            score_df = pd.read_csv(home_dir + '/' + file, sep='\t')
            score_df = score_df.drop(score_df.columns[[0]], axis = 1)
            top_ten = score_df[:10]
            top_ten["Organisms Compared to Synbio"] = top_ten["Organisms Compared to Synbio"].str.replace("_matches.tsv", "", regex=False)
            top_ten = top_ten.rename(columns={'Organisms Compared to Synbio': 'Organism Bin ID'})
            top_ten_master_list.append(top_ten)
            # print(top_ten)
            # final_hist = pd.concat(top_ten_master_list, axis=1, ignore_index=True)
            # st.dataframe(final_hist)
            st.subheader(biome_name)
            fig, ax = plt.subplots()
            norm = Normalize(vmin=top_ten["Difference Score"].min(), vmax=top_ten["Difference Score"].max())
            cmap = plt.cm.Greys_r  # reversed Greys: darker = lower

            # Map each value to a color
            colors = [cmap(norm(val)) for val in top_ten["Difference Score"]]
            ax = sns.barplot(data=top_ten, x = 'Organism Bin ID', y = 'Difference Score', ax=ax, palette='Blues_r')
            ax.bar_label(ax.containers[0], fontsize=10);
            plt.xticks(rotation=75)  

            plt.tight_layout()  
            plt.show()
            st.pyplot(fig)
    with st.expander("## How to Understand this Output"):
        st.write('Difference Score: The amount of functional overlap between your synbio organism and each organism in the biome')
        st.write('Organism Bin ID: The Joint Genome Institute\'s categorization for specific genomes in a metagenome.')
        st.write('The genome with the lowest difference score has the most functional overlap and the greatest potential for niche competition and displacement by the foreign organism')
        st.write('This program outputs the top 10 difference scores. Look in the folder you specified for all scores')
    with st.expander('# What\'s Next'):
        st.write('New Analysis: Start the Program Over and Delete all files with the biome and organism name in it')
        st.write('This will let you do another analysis without regenerating everything')
        st.write('Take a look at generated data for in depth analysis. Additionally, head on over to JGI IMG and input the Organism Bin ID for further sequencing information')
        st.write('https://img.jgi.doe.gov/cgi-bin/mer/main.cgi?section=GenomeSearch&page=searchForm')

            # st.success('Read correct number of files')    
        # else:
            # st.write('Error in threat assessment')
    
    # final_hist = pd.concat(top_ten_master_list, axis=1, ignore_index=True)
    # st.dataframe(final_hist)
    # fig, ax = plt.subplots()
    # sns.histplot(data=final_hist, x = ), ax=ax)
    # st.pyplot(fig)






