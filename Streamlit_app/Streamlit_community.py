import streamlit as st
import base64
import os

class Streamlit_community:
    pass

    @staticmethod
    def create_download_section_from_ext_link(id, name_download):
        final_link = f'<a target="_blank" rel="noopener noreferrer" href="https://drive.google.com/u/0/uc?id='+id+'&export=download">'+name_download+'</a>'
        st.markdown(final_link, unsafe_allow_html=True)


    @staticmethod #not in use
    def get_binary_file_downloader_html(bin_file, file_label='File', name_button='download'):
        with open(bin_file, 'rb') as f:
            data = f.read()
        bin_str = base64.b64encode(data).decode()
        href = f'<a href="data:application/octet-stream;base64,{bin_str}" download="{os.path.basename(bin_file)}">{name_button} {file_label}</a>'
        return href

    @staticmethod
    def get_table_download_link(df, name_of_file="dataframe.tsv", sep='\t'):
        """Generates a link allowing the data in a given pandas dataframe to be downloaded
        in:  dataframe
        out: href string
        """
        csv = df.to_csv(index=False, sep=sep)
        b64 = base64.b64encode(csv.encode()).decode()  # some strings <-> bytes conversions necessary here
        href = f'<a href="data:file/csv;base64,{b64}">' + name_of_file + '</a>'
        return href

    @staticmethod
    @st.cache(allow_output_mutation=True)
    def get_base64_of_bin_file(bin_file):
        with open(bin_file, 'rb') as f:
            data = f.read()
        return base64.b64encode(data).decode()

    @staticmethod
    def set_png_as_page_bg(png_file):
        bin_str = Streamlit_community.get_base64_of_bin_file(png_file)
        page_bg_img = '''
        <style>
        body {
        background-image: url("data:image/png;base64,%s");
        background-size: cover;
        }
        </style>
        ''' % bin_str
        st.markdown(page_bg_img, unsafe_allow_html=True)
        return

    @staticmethod
    def display_picture_with_hyperlink(picture_path,link): #figure out how to set the size of the picture
        with open(picture_path, "rb") as img_file:
            my_string = base64.b64encode(img_file.read())
        final_link = '\\\\'+link
        html = f"<a href='{final_link}'><img src='data:image/png;base64,{my_string.decode('utf-8')}'></a>"
        st.markdown(html, unsafe_allow_html=True)

    @staticmethod
    def rerun_script_from_top():
        raise st.script_runner.RerunException(st.script_request_queue.RerunData(None))