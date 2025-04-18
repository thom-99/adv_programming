import streamlit as st 
from PIL import Image 
from pathlib import Path
import faster 


image = Image.open('images/mitochondria-1.jpg')
st.image(image)


st.title("""Mitochondrial DNA Analyzer""")
st.markdown("tool developed by **Thom**, **Ale**, **Matteo** and **Fabio**")
st.markdown("the [source code](https://github.com/thom-99/adv_programming/tree/main) of the project is available on github")



# File uploader 
uploaded_file = st.file_uploader("Choose a file", type=["txt", "fasta", "fa"])

if uploaded_file is not None:
    temp_path = Path(f"temp_{uploaded_file.name}")  
    
    # Write the uploaded file to a temporary file, wb means that we write in binary mode
    # getbuffer gets all the raw byte data from the uploaded file 
    with open(temp_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    
    file = faster.Fasta(str(temp_path))
    
    fastafile = file.parse()
    st.success('file parsed successfully')


    st.markdown("## analysis section")

    #select an header, show sequence statistics
    dictionary = file.process()
    headers = dictionary.keys()
    selected_header = st.selectbox('select a sequence to analyze', options=headers)
    if selected_header:
        sequence = dictionary[selected_header]
        st.markdown(f"header: {selected_header}")
        st.markdown(f"sequence length: {len(sequence)}")
        
        analyzer = faster.Analyzer(sequence)
        gc_content = analyzer.GC_content(sequence)
        st.markdown(f"the GC content is {gc_content}")


        #alignment with another sequence
        if st.checkbox('perform alignment'):
            sequence2 = st.text_input('align the sequence with a new given one')
            if sequence2:
                st.text('below the local alignment of the two sequences is shown')
                st.text(sequence)
                aligned_seq1, aligned_seq2, score = analyzer.align(sequence, sequence2)
                match_line= ""
                for a, b in zip(aligned_seq1, aligned_seq2):
                    if a == b:
                        match_line += "|"
                    else:
                        match_line += " "
                alignment_block = f"{aligned_seq1}\n{match_line}\n{aligned_seq2}"
                st.code(alignment_block, language="text")
                st.markdown(f"**Alignment score:** `{score}`")


        #search genomic motifs in the sequence 
        if st.checkbox('search for genomic motifs'):
            motifs = st.text_input('type the genomic motifs separated by a whitespace')
            motifs_list = motifs.split(' ')
            #creating a list of sequence motif objects
            motifs_obj = [faster.Sequence_motif(x) for x in motifs_list]
            #class method call
            df = analyzer.motifs_analysis(*motifs_obj)
            st.text('motifs search results')
            st.dataframe(df)

