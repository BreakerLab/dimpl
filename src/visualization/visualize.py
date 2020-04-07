import plotly.graph_objs as go
from src.models.train_model import genome_svm_selection
from IPython.display import display

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

