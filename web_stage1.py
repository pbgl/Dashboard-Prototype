import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import pandas as pd
import numpy as np
import dash_table
from dash.dependencies import Input, Output
import dash_daq as daq

'''my_file=open('format.css','r')
external_stylesheets = my_file.readlines()'''
x={'label':'all','value':'all'}
tabs_styles = {
    'height': '85px'
}
tab_style = {
    'borderBottom': '1px solid #d6d6d6',
    'padding': '6px',
    'fontWeight': 'bold'
}

tab_selected_style = {
    'borderTop': '1px solid #d6d6d6',
    'borderBottom': '1px solid #d6d6d6',
    'backgroundColor': '#119DFF',
    'color': 'white',
    'padding': '6px'
}
app = dash.Dash(__name__)

# Importing data
SnpSiftData=pd.read_csv(r'data/snpsiftdata.csv',delimiter='\t',encoding='UTF-8')
chromosome_name=pd.read_csv(r'data/chromosome_name_mapping.csv',delimiter='\t',encoding='UTF-8')
Genotype_Data=pd.read_csv(r'data/Genotype_Data.csv',delimiter='\t',encoding='UTF-8')
passport=pd.read_csv(r'data/passport.csv',delimiter='\t',encoding='UTF-8')

passport.replace(np.NaN,'NA')

Genotype_Data_1=pd.melt(Genotype_Data,id_vars=['CHROM','POS'],var_name='Sample_ID',value_name='GT')
Genotype_Data_1['Chrom_Pos']=Genotype_Data_1['CHROM'].astype(str).str.cat(Genotype_Data_1['POS'].astype(str),sep='_')
Genotype_Data_1.insert(4,'Variety',Genotype_Data_1['Sample_ID'].map(passport.set_index('Sample-ID')['Variety']))
Genotype_Data_1.insert(5,'Generation',Genotype_Data_1['Sample_ID'].map(passport.set_index('Sample-ID')['Generation']))  
SnpSiftData['Chrom_Pos']=SnpSiftData['CHROM'].astype(str).str.cat(SnpSiftData['POS'].astype(str),sep='_')
SnpSiftData.insert(5,'chrome_name',SnpSiftData['CHROM'].map(chromosome_name.set_index('Contig')['Chromosome']))
clicker=0
variants=['complex','snp','mnp','del','ins']
impact=['HIGH','LOW','MODERATE','MODIFIER']
# assume you have a "wide-form" data frame with no index
# see https://plotly.com/python/wide-form/ for more options
df = pd.DataFrame({"x": [1, 2, 3], "SF": [4, 1, 2], "Montreal": [2, 4, 5]})

#fig = px.bar(df, x="x", y=["SF", "Montreal"], barmode="group")
variant_options=[{'label':i, 'value':i} for i in SnpSiftData.TYPE.unique() ]
variant_options.insert(0,x)
impact_options=[{'label':i, 'value':i} for i in SnpSiftData['ANN[*].IMPACT'].unique()]
impact_options.insert(0,x)
effect_options=[{'label':i, 'value':i} for i in SnpSiftData['ANN[*].EFFECT'].unique()]
effect_options.insert(0,x)

app.layout =html.Div(children=[
    html.Div(className='container',children=[
        
		html.Div(className='left_container',children=[        
			html.Div(children=[
				html.H2(className='h2',children='FAO/IAEA-PBGL',style={'color': '#056aae'}),
                html.H3(className='h3',children='Coffee Mutants Browser',style={'color': '#056aae'})
			]),
			html.Div(children=[
				html.Div(children=[
					html.Div(children=[
						html.Div(children=[
							html.H5("Variant Search Parameters")
									#style={'width': '40%', 'float': 'left'})
							])
						])
					]),
				html.Div(children=[
					html.Div(children=[
						dcc.Tabs(id='tabs_id',value='Gene_tab',children=[
							dcc.Tab(label='Gene Identifier',value='Gene_tab',children=[
								html.Div(children=[
									html.H5(children='Gene Name:',style={'color':'red'}),
									dcc.Input(id='Gene_Identifier',placeholder='Type to pick some TAIR IDs...',value='all',style={'marginBottom': '1.5em'})
                                    ]),
                                html.Div(children=[
                                    html.H5(children='Max Distance from Gene (bp):',style={'color':'red'}),
                                    dcc.Input(id='Distance',placeholder='Max Distance from the gene',value='',
                                              style={'marginBottom': '1.5em'},type='number')
                                ])
								
							],style=tab_style, selected_style=tab_selected_style),
							dcc.Tab(label='Range',value='range_tab',children=[
								html.Div(children=[
									html.H5(children='Chromosome Name:',style={'color':'red'}),
									dcc.Input(id='chromosome_name',type='text',placeholder='Chromosome Name',
											  value='all')
									]),
								html.Div(children=[
									html.H5( children='Start:',style={'color':'red'}),
									dcc.Input(id='start_pos',type='number',placeholder='Start Position',
											  value='all')
									]),
								html.Div(children=[
									html.H5(children='End:',style={'color':'red'}),
									dcc.Input(id='end_pos',type='number',placeholder='End Position',
											  value='all',style={'marginBottom': '1.5em'})
									]),
								],style=tab_style, selected_style=tab_selected_style),
							dcc.Tab(label='Position',value='Pos_tab',children=[
								html.Div(children=[
									html.H5(children='Chromosome Name:',style={'color':'red'}),
									dcc.Input(id='chromosome_name_pos',type='text',placeholder='Chromosome Name')
									]),
								html.Div(children=[
									html.H5(children='Position',style={'color':'red'}),
									dcc.Input(id='position',type='number',placeholder=' Position',style={'marginBottom': '1.5em'})
									]),
								],style=tab_style, selected_style=tab_selected_style)
						])
					])
				]),#,style={'width': '20%', 'float': 'left'}),
                
				html.Div(children=[
                    html.H5('Variant Filter'),
                    dcc.Tabs(value='variants_tab',children=[
                        dcc.Tab(label='Variant Type',value='variants_tab',children=[
                            dcc.Checklist(id='variant',
                                          options=[{'label':i,'value':i} for i in variants],
                                          value=['complex','snp','mnp','del','ins'])
                        ],style=tab_style, selected_style=tab_selected_style),
                        dcc.Tab(label='Impact Type',value='impact_tab',children=[
                            dcc.Checklist(id='impact',
                                          options=[{'label':i,'value':i} for i in impact],
                                          value=['HIGH','LOW','MODERATE','MODIFIER'])
                        ],style=tab_style, selected_style=tab_selected_style),
                        dcc.Tab(label='Effect Type',value='effect_tab',children=[
                            dcc.Dropdown(id='effect_type',
                                         options=effect_options,
                                         value='all')
                        ],style=tab_style, selected_style=tab_selected_style),
                    ])
				],style={"border":"3px red dotted"}),
                html.Div(children=[
                    html.H5('Plant Filter'),
                    dcc.Tabs(value='plant_tab',children=[
                        dcc.Tab(label='Variety',value='variety_tab',children=[
                            dcc.Checklist(id='variety',
                                          options=[{'label':i,'value':i} for i in Genotype_Data_1.Variety.unique()],
                                          value=['Venetia','Caturra','Catuai'])
                        ],style=tab_style, selected_style=tab_selected_style),
                        dcc.Tab(label='Generation',value='generation_tab',children=[
                            dcc.Checklist(id='generation',
                                         options=[{'label':i,'value':i} for i in Genotype_Data_1.Generation.unique()],
                                         value=['M0'])
                        ],style=tab_style, selected_style=tab_selected_style),
                    ])
				],style={"border":"3px blue dotted"}),
                html.Div(className='control-tab',children=[
                    html.H5('Noise Removal'),
                    html.Div(children=[
                        html.Div(children='Plant With REF Allel (0/0)'),
                        daq.ToggleSwitch(
                                    id='noisy_snps',
                                    label=['hide', 'show'],
                                    color='#009DFF',
                                    size=35,
                                    value=False
                                ),
                    ],style={'inline':True}),
                    html.Div(children=[
                        html.Div(children='Multi Allelic'),
                        daq.ToggleSwitch(
                                    id='Multi_allelic',
                                    label=['hide', 'show'],
                                    color='#009DFF',
                                    size=35,
                                    value=True
                                ),
                    ],style={'inline':True})
                ],style={"border":"3px black dotted"}),
                
                html.Div(className='right',children=[
                    html.Button(id='search_button',n_clicks=0, children='Search',style={'backgroundColor':'#FF0000','color':'black','fontWeight': 'bold',
                                                                                        'padding-left':'25%', 'padding-right':'25%','float':'center'})
                ]),
            ])
        ]),
        html.Div(className='right_container',children=[
            html.Div(children=[
                html.Div(children=[
                    dash_table.DataTable(id='table_data',
                                         columns=[
                                             {'name':'Gene','id':'ANN[*].GENE','type':'text'},
                                             {'name':'Chromosome','id':'chrome_name','type':'text'},
                                             {'name':'Contig Name','id':'CHROM_x','type':'text'},
                                             {'name':'Position','id':'POS_y','type':'text'},
                                             {'name':'Sample_ID','id':'Sample_ID','type':'text'},
                                             {'name':'Variety','id':'Variety','type':'text'},
                                             {'name':'Generation','id':'Generation','type':'text'},
                                             {'name':'GT','id':'GT','type':'text'},
                                             {'name':'REF','id':'REF','type':'text'},
                                             {'name':'ALT','id':'ALT','type':'text'},
                                             {'name':'TYPE','id':'TYPE','type':'text'},
                                             {'name':'IMPACT','id':'ANN[*].IMPACT','type':'text'},
                                             {'name':'EFFECT','id':'ANN[*].EFFECT','type':'text'},
                                             {'name':'Distance','id':'ANN[*].DISTANCE','type':'text'},
                                             {'name':'ID','id':'ID','type':'text'},
                                             ],
                                         sort_action="native",
                                         sort_mode="multi",
                                         page_action='none',
                                         fixed_rows={'headers': True,'data':0},
                                         style_table={'maxheight': '1500','overflowY':'scroll'},
                                         filter_action='native',
                                         style_data={'width': '{}%'.format(100. / 15), # 14 is the number of columns to display on the page as a result.
                                                     'border': '1px solid #0000FF' },
                                         style_header={
                                             'backgroundColor': '#C0C0C0',
                                             'fontcolor':'#0000FF',
                                             'fontWeight': 'bold'},
                                         style_cell_conditional=[{'textAlign': 'left'}],
                                         )
                    ],style={'height':'100%'})
                ])
            ])
        ]),
    ])

def filters(test_data,variant_value,impact_value,effect_value,Multi_allelic_value):
    test_data_1=pd.DataFrame()
    test_data_2=pd.DataFrame()
    test_data_3=pd.DataFrame()
    if Multi_allelic_value == True:
        for i in variant_value:
            test_data_1=test_data_1.append(test_data[test_data['TYPE'].str.contains(i)])
    
    if Multi_allelic_value == False:
        for i in variant_value:
            test_data_1=test_data_1.append(test_data[test_data['TYPE']==i])
        
        
    
    if len(impact_value) == 4:
        test_data_2=test_data_1
        
    else:
        for j in impact_value:
            test_data_2=test_data_2.append(test_data_1[test_data_1['ANN[*].IMPACT']==j])
            print(test_data_2)
        
    if effect_value == 'all':
        test_data_3=test_data_2
        print(test_data_3.head())
    else:
        test_data_3=test_data_2[test_data_2['ANN[*].EFFECT']==effect_value]
        print(test_data_3.head())
        
    return test_data_3


@app.callback(
    Output('table_data',"data"),
    [Input('tabs_id',"value"),
     Input('Gene_Identifier',"value"),
     Input('chromosome_name',"value"),
     Input('start_pos',"value"),
     Input('end_pos',"value"),
     Input('chromosome_name_pos',"value"),
     Input('position',"value"),
     Input('variant',"value"),
     Input('impact',"value"),
     Input('effect_type',"value"),
     Input('search_button',"n_clicks"),
     Input('noisy_snps',"value"),
     Input('Multi_allelic',"value"),
     Input('variety',"value"),
     Input('generation',"value"),
     Input('Distance',"value"),])
def GenoFiltering(tabs_value,Geno_value,chrome_name_value,start_pos_value,end_pos_value,choromosome_name_third_tab,pos_third_tab,
                  variant_value,impact_value,effect_value,n_clicks,noisy_value,Multi_allelic_value,variety_value,generation_value,distance_value):
    global clicker
    
    final_data_3=pd.DataFrame()    
    if n_clicks > clicker:
        if tabs_value == 'Gene_tab':
            #if Geno_value == 'all' and chrome_name_value=='all' and start_pos_value == 'all' and end_pos_value == 'all':
            #test_data=SnpSiftData
            #filtered_data=filters(test_data,variant_value,impact_value,effect_value)
            test_data=SnpSiftData[SnpSiftData['ANN[*].GENE']==Geno_value]
            if distance_value == '':
                test_data=test_data
            else:
                test_data=test_data[test_data['ANN[*].DISTANCE'] <= distance_value]
                
            filtered_data=filters(test_data,variant_value,impact_value,effect_value,Multi_allelic_value)
            print('Gene Tab Initiated')
        if tabs_value == 'range_tab':
            
            test_data=SnpSiftData[SnpSiftData['chrome_name']==chrome_name_value]
            if start_pos_value != 'all' and end_pos_value != 'all':
                test_data_1=test_data[tart_pos_value < test_data['POS'] and test_data['POS'] < end_pos_value]
            else:
                test_data_1=test_data
            filtered_data=filters(test_data_1,variant_value,impact_value,effect_value,Multi_allelic_value)
                
        if tabs_value == 'Pos_tab':
            test_data=SnpSiftData[SnpSiftData['chrome_name']==choromosome_name_third_tab]
            test_data_1=test_data[test_data['POS'] == pos_third_tab]
            filtered_data=filters(test_data_1,variant_value,impact_value,effect_value,Multi_allelic_value)
        
        final_data=pd.merge(Genotype_Data_1,filtered_data,on='Chrom_Pos')
        final_data_1=final_data[['ANN[*].GENE','chrome_name','CHROM_x','POS_y','Sample_ID','Variety','Generation','GT','REF','ALT','TYPE','ANN[*].IMPACT','ANN[*].EFFECT','ANN[*].DISTANCE','ID']]
        clicker=clicker+1
        print(final_data_1.head(10))
        if noisy_value == True:
            final_data_2=final_data_1
        else:
            print(final_data_1.head(10))
            final_data_2=final_data_1[final_data_1['GT'] != '0/0']
            
        for i in variety_value:
            final_data_3=final_data_3.append(final_data_2[final_data_2['Variety']==i])
        
            
    else:
        final_data_3=pd.DataFrame(columns=[['ANN[*].GENE','chrome_name','CHROM_x','POS_y','Sample_ID','Variety','Generation','GT','REF','ALT','TYPE','ANN[*].IMPACT','ANN[*].EFFECT','ANN[*].DISTANCE','ID']],data=None)
    #final_data_3.replace(np.NaN,'NA')
    return final_data_3.to_dict('records')


    
    
    

if __name__ == '__main__':
    app.run_server(debug=True)
