import sys
import os
import pandas as pd
import plotly.graph_objs as go
import tarfile
import pickle
import plotly as py
import shutil
from sklearn import svm
from ipywidgets import interactive
from src.models.train_model import genome_svm_selection
from IPython.display import display
from src.data.rfam_db import rfam_session, Genome
from src.data.make_dataset import extract_igrs, annotate_igrs, download_genome
from Bio import SeqIO, SeqRecord

def graph_layout(genome):

    ytickvals = list(range(0, 110, 10))
    yticktext = ['<b>0</b>', '', '<b>20</b>', '', '<b>40</b>', '', '<b>60</b>', '', '<b>80</b>', '', '<b>100</b>']
    xtickvals = list(range(10, 100, 10)) + list(range(100, 1000, 100)) + list(range(1000, 11000, 1000))
    xticktext = ['<b>10</b>'] + [''] * 8 + ['<b>100</b>'] + [''] * 8 + ['<b>1000</b>'] + [''] * 8 + ['<b>10000</b>']

    layout = go.Layout(title=dict(text="<b><i>{}</i></b>".format(genome.scientific_name), x=0.45, y=0.925),
                       font=dict(family="Arial", size=18, color='black'),
                       yaxis=dict(title="<b> GC Content (%) </b>",
                                  showline=True, showgrid=False,
                                  linecolor='black',
                                  linewidth=2, tickwidth=2,
                                  tickfont=dict(size=16),
                                  ticks='outside',
                                  tickvals=ytickvals, ticktext=yticktext,
                                  mirror=True,
                                  range=[0, 100]),
                       xaxis=dict(title=dict(text="<b> IGR Length </b>"),
                                  type="log",
                                  tickangle=0,
                                  tickvals=xtickvals, ticktext=xticktext,
                                  showline=True, showgrid=False,
                                  linecolor='black', linewidth=2,
                                  tickfont=dict(size=16), tickwidth=2,
                                  ticks='outside',
                                  mirror=True,
                                  range=[1, 4],),
                       width=900,
                       hovermode='closest',
                       legend=dict(y=0.5, x=1.05, borderwidth=2),
                       height=600,
                       autosize=False,
                       margin=dict(autoexpand=False,l=75,r=250,b=100,t=100,pad=5),
                       plot_bgcolor="white",paper_bgcolor='white')

    return layout

def graph_genome(annotated_df, selection=None):

    annotated_df = annotated_df.copy()

    category_style = {
        "No Known RNA": ["triangle-up-open", "rgba(192,192,192,0.9)", 'skip'],
        "Selected IGR": ["triangle-up", "rgba(192,192,192,0.9)", 'skip'],
        "sRNA": ["triangle-down", "rgba(192,192,192,0.9)", 'text'],
        "tRNA": ["diamond", "rgba(55,126,184,0.8)", 'text'],
        "rRNA": ["square", "rgba(55,126,184, 0.8)", 'text'],
        "Intron": ["cross", "rgba(255,255,51, 0.8)", 'text'],
        "RNase P": ["triangle-right", "rgba(255,127,0, 0.8)", 'text'],
        "6S RNA": ["triangle-left", "rgba(55,126,184,0.8)", 'text'],
        "tmRNA": ["diamond", "rgba(255,127,0,0.8)", 'text'],
        "Riboswitch": ["star", "rgba(228,26,28,0.8)", 'text'],
        "Ribozyme": ["x", "rgba(247,129,191,0.8)", 'text'],
        "Miscellaneous": ["pentagon", "rgba(166,86,40,0.8)", 'text'],
        "Multiple": ["star-diamond", "rgba(152,78,163,0.8)", 'text']
    }


    knowns = annotated_df['category'] != 'No Known RNA'
    total_igrs = len(knowns)
    total_knowns = sum(knowns)


    if selection is not None:

        # Set the values of category to "Selected IGR" vs
        annotated_df.loc[annotated_df["rfam_acc"].isnull() & selection, "category"] = "Selected IGR"
        annotated_df.loc[annotated_df["rfam_acc"].isnull() & ~selection, "category"] = "No Known RNA"


        unique_categories = annotated_df["category"].unique()


        knowns_included = sum(knowns & selection)
        unknowns_included = sum(~knowns & selection)
        print("Number of known IGRs included:   {} ({:.1%})".format(knowns_included,
                                                                            knowns_included / total_knowns))
        print("Number of unknown IGRs included: {} ({:.1%})".format(unknowns_included, unknowns_included/total_igrs))
        print("Fold Enrichment: {:5.2f}".format((knowns_included/sum(selection))/(total_knowns/total_igrs) ))

        annotated_df['selection'] = selection
        point_selection = [
            list(annotated_df[annotated_df['category'] == category].reset_index().query('selection').index) for
            category in unique_categories]
    else:

        unique_categories = annotated_df["category"].unique()
        point_selection = [None] * len(unique_categories)

    scatter_plots = [go.Scatter(y=annotated_df[annotated_df['category'] == category]['gc'],
                                x=annotated_df[annotated_df['category'] == category]['length'],
                                name=category,
                                mode='markers',
                                selectedpoints= point_selection[index],
                                text=annotated_df[annotated_df['category'] == category]['description'],
                                hoverinfo=category_style[category][2],
                                marker=dict(size=10,
                                            symbol=category_style[category][0],
                                            color=category_style[category][1])
                                ) for index, category in enumerate(unique_categories)]

    return scatter_plots


def interactive_selection(annotated_df, layout, gamma=0.001, class_weight_mod=1):

    selection = genome_svm_selection(annotated_df, gamma=gamma, class_weight_mod=class_weight_mod)

    # Build the plotly scatter plots for each category
    scatter_plots = graph_genome(annotated_df, selection=selection)

    fig = go.FigureWidget(data=scatter_plots, layout=layout)

    display(fig)

    return fig

def display_genome(upid):
    
    session = rfam_session()
    genome =  session.query(Genome).get(upid)
    session.close()
    download_genome(genome)
    igr_df = extract_igrs(genome, igr_length_cutoff=1)
    annotated_df = annotate_igrs(genome, igr_df)
    scatter_plots = graph_genome(annotated_df)
    layout = graph_layout(genome)
    fig = go.FigureWidget(data=scatter_plots, layout=layout)
 
    return annotated_df, fig, layout, genome

def prepare_selection(annotated_df):
    y = (annotated_df['category'] != 'No Known RNA') & (annotated_df['category'] != 'sRNA') 
    total_igrs = len(y)
    total_knowns = y.sum()
    total_unknowns = total_igrs - total_knowns 
    
    return (y, total_igrs, total_knowns, total_unknowns)

def build_interactive_fn (annotated_df, layout, genome):
    
    y, total_igrs, total_knowns, total_unknowns = prepare_selection(annotated_df)    

    def interactive_fn(class_weight_mod=0.5,  c_exp=2, gamma_exp=-2,):
        class_weight = {False: total_knowns / total_igrs, True: (total_unknowns / total_igrs * class_weight_mod)}
        svm_clf = svm.SVC(C=10**c_exp, class_weight=class_weight, gamma=10**(gamma_exp), random_state=0)
        svm_clf.fit(annotated_df.loc[:, ["gc", "log_length"]], y)
        selection = pd.Series(svm_clf.predict(annotated_df.loc[:, ["gc", "log_length"]]))
        scatter_plots = graph_genome(annotated_df, selection=selection)
        fig = go.FigureWidget(data=scatter_plots, layout=layout)
        display(fig)
    
    return interactive_fn

def save_selected_IGRs(interactive_plot, annotated_df, genome):
    class_weight_mod = interactive_plot.kwargs["class_weight_mod"]
    c_exp = interactive_plot.kwargs["c_exp"]
    gamma_exp = interactive_plot.kwargs["gamma_exp"]

    output_folder="data/interim/{}/selection_{}_{}_{}".format(genome.assembly_acc, class_weight_mod, c_exp, gamma_exp)
    if not os.path.exists(output_folder + '/data_files'):
        os.makedirs(output_folder + '/data_files')

    # Re-create the selection
    y, total_igrs, total_knowns, total_unknowns = prepare_selection(annotated_df)    
    class_weight = {False: total_knowns / total_igrs, True: (total_unknowns / total_igrs * class_weight_mod)}
    svm_clf = svm.SVC(C=10**c_exp, class_weight=class_weight, gamma=10**gamma_exp, probability=True, random_state=0)
    svm_clf.fit(annotated_df.loc[:, ["gc", "log_length"]], y)

    # Save the selection classifier to a pickle                                                          
    svm_pickle = pickle.dumps(svm_clf)
    with open("{}/data_files/svmclf.pickle".format(output_folder,genome.assembly_acc, class_weight_mod, c_exp, gamma_exp), 'wb') as svm_pickle_file:
        svm_pickle_file.write(svm_pickle)
    selection = pd.Series(svm_clf.predict(annotated_df.loc[:, ["gc", "log_length"]]))

    # Save a graph of the genome.
    scatter_plots = graph_genome(annotated_df, selection=selection)
    layout = graph_layout(genome)
    fig = go.FigureWidget(data=scatter_plots, layout=layout)
    fig.write_image("{}/data_files/{}_plot.svg".format(output_folder,genome.scientific_name.replace(' ','_')))
    py.io.write_json(fig, "{}/data_files/{}_json.plotly".format(output_folder,genome.scientific_name.replace(' ','_')))

    selected_unknowns = selection & (annotated_df['category'] == 'No Known RNA')

    # Save a fasta file with all the selected IGRs
    selected_igr_list = [SeqRecord.SeqRecord(row.sequence, id=("{}/{}-{}".format(row.accession, row.start +1, row.end))) 
                         for row in annotated_df.loc[selected_unknowns, ["accession","start","end","sequence"]].itertuples()]

    if not os.path.exists(output_folder + '/igr_fastas'):
        os.makedirs(output_folder + '/igr_fastas')

    for igr in selected_igr_list:
        outputfilename = "{}/igr_fastas/{}.fasta".format(output_folder, igr.id.replace('/','_'))
        SeqIO.write(igr, outputfilename, "fasta")

    annotated_df.to_csv("{}/data_files/annotated_igrs.csv".format(output_folder), index=False)
    
    #Block 6
    
    if not os.path.exists(output_folder + '/scripts'):
        os.makedirs(output_folder + '/scripts')

    shutil.copy('src/shell/cluster.conf', '{}/scripts'.format(output_folder))
    shutil.copy('src/shell/make_tar.sh', '{}/scripts'.format(output_folder))
    shutil.copy('src/shell/blast_source_template.sh', '{}/scripts/blast_source.sh'.format(output_folder))
    shutil.copy('src/shell/blast_run_template.sh', '{}/blast_run.sh'.format(output_folder))

    if not os.path.exists(output_folder + '/blast_xml'):
        os.makedirs(output_folder + '/blast_xml')
    if not os.path.exists(output_folder + '/output'):
        os.makedirs(output_folder + '/output')

    with open("{}/scripts/blast_jobfile.sh".format(output_folder), 'w') as jobfile:
        for igr in selected_igr_list:
            fasta_filename = "igr_fastas/{}.fasta".format(igr.id.replace('/','_'))
            xml_filename =  "blast_xml/{}.xml".format(igr.id.replace('/','_'))
            jobfile.write("source scripts/blast_source.sh; $BLASTCMD {} > {}\n".format(fasta_filename, xml_filename))

    with tarfile.open("data/export/{}_{}_selection_{}_{}_{}_blastdata.tar.gz".format('_'.join(genome.scientific_name.split(' ')[0:2]), genome.assembly_acc, class_weight_mod, c_exp, gamma_exp), "w:gz") as tar:
        tar.add(output_folder, arcname="{}_{}_selection_{}_{}_{}_blastdata".format('_'.join(genome.scientific_name.split(' ')[0:2]), genome.assembly_acc, class_weight_mod, c_exp, gamma_exp))
        print("\nTarfile created:",tar.name)
    return
