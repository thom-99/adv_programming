# Faster
a web tool for processing and analyzing fasta files 

The front-end is developed integrally in Python using [streamlit](https://streamlit.io/). 

To install streamlit run one of the following commands in your terminal:
```
pip install streamlit
conda install streamlit   #if you prefer to use a conda enviroment 
```
To run the app navigate to the folder containing the app and run the following command:
```
streamlit run app.py
```

## DEMO

1. after running the app, you will be asked to input a file from your machine. It is recommended for your own sake that you store your fasta file using a .fasta or .fa extension, however the app also supports the upload of text files.

  ![Screenshot 2025-04-15 150032](https://github.com/user-attachments/assets/b84c8457-9791-4cb6-802b-42045a0c3f8d)

2. then you need only to navigate to where your file of interest is, you can use the fasta file in the /example_data folder as a demo
   
  ![Screenshot 2025-04-15 150112](https://github.com/user-attachments/assets/3601bd1b-e5ad-4298-9c6a-ba3d1b057a40)

3. if there are no errors a green banner will appear, and you will be able to proceed with the analysis of your file. In this section you have to select a header, which corresponds to a sequence of your file, and after that you can use advanced features such as motif search or alignment.
   
   ![Screenshot 2025-04-15 150157](https://github.com/user-attachments/assets/8f09fbf7-aabd-4e34-a7d3-ae26d2e935a2)

