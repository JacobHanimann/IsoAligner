import streamlit as st
import base64
import os
import uuid
import re

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
    def download_button(object_to_download, download_filename, button_text,df, sep):
        """
        Generates a link to download the given object_to_download.

        From: https://discuss.streamlit.io/t/a-download-button-with-custom-css/4220
        Params:
        ------
        object_to_download:  The object to be downloaded.
        download_filename (str): filename and extension of file. e.g. mydata.csv,
        some_txt_output.txt download_link_text (str): Text to display for download
        link.
        button_text (str): Text to display on download button (e.g. 'click here to download file')
        pickle_it (bool): If True, pickle file.
        Returns:
        -------
        (str): the anchor tag to download object_to_download
        Examples:
        --------
        download_link(your_df, 'YOUR_DF.csv', 'Click to download data!')
        download_link(your_str, 'YOUR_STRING.txt', 'Click to download text!')
        """
        # if pickle_it:
        #     try:
        #         object_to_download = pickle.dumps(object_to_download)
        #     except pickle.PicklingError as e:
        #         st.write(e)
        #         return None

        # else:
        #     if isinstance(object_to_download, bytes):
        #         pass

        #     elif isinstance(object_to_download, pd.DataFrame):
        #         object_to_download = object_to_download.to_csv(index=False)

        #     # Try JSON encode for everything else
        #     else:
        #         object_to_download = json.dumps(object_to_download)

        try:
            # some strings <-> bytes conversions necessary here
            csv = df.to_csv(index=False, sep=sep)
            b64 = base64.b64encode(csv.encode()).decode()  # some strings <-> bytes conversions necessary here
        except AttributeError as e:
            b64 = base64.b64encode(object_to_download).decode()

        button_uuid = str(uuid.uuid4()).replace("-", "")
        button_id = re.sub("\d+", "", button_uuid)

        custom_css = f""" 
            <style>
                #{button_id} {{
                    display: inline-flex;
                    align-items: center;
                    justify-content: center;
                    background-color: rgb(255, 255, 255);
                    color: rgb(38, 39, 48);
                    padding: .25rem .75rem;
                    position: relative;
                    text-decoration: none;
                    border-radius: 4px;
                    border-width: 1px;
                    border-style: solid;
                    border-color: rgb(230, 234, 241);
                    border-image: initial;
                }} 
                #{button_id}:hover {{
                    border-color: rgb(246, 51, 102);
                    color: rgb(246, 51, 102);
                }}
                #{button_id}:active {{
                    box-shadow: none;
                    background-color: rgb(246, 51, 102);
                    color: white;
                    }}
            </style> """

        dl_link = (
                custom_css
                + f'<a download="{download_filename}" id="{button_id}" href="data:file/txt;base64,{b64}">{button_text}</a><br><br>'
        )
        # dl_link = f'<a download="{download_filename}" id="{button_id}" href="data:file/txt;base64,{b64}"><input type="button" kind="primary" value="{button_text}"></a><br></br>'

        st.markdown(dl_link, unsafe_allow_html=True)

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