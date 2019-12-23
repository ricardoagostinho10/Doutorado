from flask import Flask, render_template, request, url_for, flash, redirect, session, jsonify
import numpy as np
from cardiojunction import cardiojunction
import warnings
warnings.filterwarnings('ignore')
app = Flask(__name__)
app.secret_key = 'fsdfsdfsdfsdfsfsd'
import plotly.plotly as py
import plotly.graph_objs as go
import pandas as pd
import json
import plotly

application = app

@app.route('/teste')
def get_plot():

    protocol = int(request.args.get('selProtocolo'));
        
    Model_ICaL = int(request.args.get('selCalcio'));

    Model_Na = int(request.args.get('selSodio'));

    Model_Ito =  int(request.args.get('selPotassioTrans'));

    Model_IKr = int(request.args.get('selPotassioR'));

    Model_IKs =  int(request.args.get('selPotassioS'));

    Model_Force = int(request.args.get('selForca'));

    cellLength  = int(request.args.get('txtComCel'));

    Lsarc = float(request.args.get('txtComSar'));
	
    BLOCKSRPUMP = float(request.args.get('txtBloqueioPump'))/100;
    STIMULUSRPUMP = float(request.args.get('txtEstimuloPump'))/100;
    BLOCKNCX = float(request.args.get('txtBloqueioNCX'))/100;
    STIMULUSNCX = float(request.args.get('txtEstimuloNCX'))/100;
    CAFEINA = int(request.args.get('selCafeina'))/100;
    BLOCKCICR = float(request.args.get('txtBloqueioCICR'))/100;
    BLOCKIKs = float(request.args.get('txtBloqueioIKs'))/100;
    BLOCKIKr = float(request.args.get('txtBloqueioIKr'))/100;
    BLOCKItof = float(request.args.get('txtBloqueioIKtof'))/100;
    BLOCKItos = float(request.args.get('txtBloqueioIKtos'))/100;
    BLOCKINa = float(request.args.get('txtBloqueioINaFast'))/100;
    BLOCKICaL = float(request.args.get('txtBloqueioICaL'))/100;
    STIMULUSICaL = float(request.args.get('txtEstimuladorICaL'))/100;
   
    t_ap = int(request.args.get('txtInstanteAplicacaoPulso'));
       
    L = int(request.args.get('txtTempoSimulacao'));

    L = L- t_ap;

    Ap = int(request.args.get('txtEstimuloVoltagem'));
        
    v_resting = int(request.args.get('txtPotencialRepouso'));

    tap = int(request.args.get('txtTempoSimulacaoAplicacao'));

    w = int(request.args.get('txtLarguraPulso'));

    f = int(request.args.get('txtFrequencia'));
        
    Delay = 1000*1/f - w;

    A_inj = float(request.args.get('txtAmplitudeInjetada'));
        
    tig = int(request.args.get('txtInstantePlotagem'));
	
    c = cardiojunction();
	
    [t,I_app,v,I_Ca,I_Na,z37,z35,z36,z30,I_kr,I_ks,C0KsC1Ks,O1KsO2Ks,OKr,C1KrC2KrC3Kr,I_kp,I_ki,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,O1,O2,J_SRCarel,J_serca,J_SRleak,z13,z15,z14,    canalSerca,I_ncx,I_nak,z34,z32,z31,z33,DifNa_sl_cl,I_cabk_junc,I_cabk_sl,I_cabk,I_pca_junc,I_pca_sl,I_pca,    Po,z40,z41,z42,z43,z44,z45,z4041,z4243,z4445,z3,z4,ONa,z484950,z515253,z5455,    z0,z1,z2,I_tos,I_tof,Of,If,Os,Is,C0fC1fC2fC3f,CI0fCI1fCI2fCI3f,C0sC1sC2sC3s,CI0sCI1sCI2sCI3s,z7,    z8,z9,z10,FORCA,Lsim,Cacy,CelL,Fcontr,N0,N1,P0,P1,P2,P3,SL,IKur,IKss]=c.principal(protocol, Model_ICaL, Model_Na, Model_Ito, Model_IKr, Model_IKs, Model_Force,cellLength, Lsarc, BLOCKSRPUMP,STIMULUSRPUMP,BLOCKNCX,STIMULUSNCX,CAFEINA,BLOCKCICR,BLOCKIKs,BLOCKIKr,BLOCKItof,BLOCKItos,BLOCKINa,BLOCKICaL,STIMULUSICaL,t_ap,L,Ap,v_resting,tap,w,f,Delay,A_inj,tig);

    #UM = np.linspace(1,1,len(t));
    #UM = UM.transpose();


    graphs = [
        dict(data=[dict(x=t,y=I_app,type='scatter')],layout=dict(title='Corrente Injetada',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2,),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=v,type='scatter'),],layout=dict(title='Voltagem',yaxis=dict(title= 'mV',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=I_Ca,type='scatter'),],layout=dict(title='ICal Current',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
		dict(data=[dict(x=t,y=I_Na,type='scatter'),],layout=dict(title='Fast Na Current',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=z37,type='scatter')],layout=dict(title='[CA]cleft',yaxis=dict(title= '[CA2+]cleft in mM',ticklen= 5,gridwidth= 2,),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=z35,type='scatter'),],layout=dict(title='[CA]ss',yaxis=dict(title= '[CA2+]ss in mM',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=z36,type='scatter'),],layout=dict(title='[CA]cyt',yaxis=dict(title= '[CA2+]cyt in mM',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
		dict(data=[dict(x=t,y=z30,type='scatter'),],layout=dict(title='[CA]SR',yaxis=dict(title= '[CA2+]SR in mM',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
		dict(data=[dict(x=t,y=I_kr,type='scatter',name='Ikr'),dict(x=t,y=I_ks,type='scatter',name='Iks')],layout=dict(title='Currents IKr and IKs',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=C0KsC1Ks,type='scatter',name='CLOSED'),dict(x=t,y=O1KsO2Ks,type='scatter',name='OPEN')],layout=dict(title='States Channel IKs',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2,),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=C1KrC2KrC3Kr,type='scatter',name='CLOSED'),dict(x=t,y=I_kr,type='scatter',name='INACTIVED'),dict(x=t,y=OKr,type='scatter',name='OPEN')],layout=dict(title='[CA]ss',yaxis=dict(title= 'States channel IKr',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=I_kp,type='scatter'),],layout=dict(title='IKp Current',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=I_ki,type='scatter'),],layout=dict(title='IKi Current',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=J_SRCarel,type='scatter')],layout=dict(title='CICR',yaxis=dict(title= 'mM/ms',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=J_serca,type='scatter')],layout=dict(title='Serca Pump SR',yaxis=dict(title= 'mM/ms',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=J_SRleak,type='scatter')],layout=dict(title='Leak SR',yaxis=dict(title= 'mM/ms',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=z13,type='scatter',name='Closed'),dict(x=t,y=1-(z14+z15+z13),type='scatter',name='Inactived - Closed')],layout=dict(title='IKp Current',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=z15,type='scatter',name='Inactived'),dict(x=t,y=z14,type='scatter',name='Open')],layout=dict(title='IKp Current',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=I_ncx,type='scatter')],layout=dict(title='NCX',yaxis=dict(title= 'mM/ms',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=I_nak,type='scatter')],layout=dict(title='NaK',yaxis=dict(title= 'mM/ms',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=z34,type='scatter')],layout=dict(title='K Concentration',yaxis=dict(title= 'mM/ms',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=z32,type='scatter',name='Na Concentration in Subsarcollema'),dict(x=t,y=z31,type='scatter',name='Na Concentration in Cleft'),dict(x=t,y=z33,type='scatter',name='Na Concentration in Cytosol')],layout=dict(title='Na Concentration',yaxis=dict(title= 'mM',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=DifNa_sl_cl,type='scatter',)],layout=dict(title='Difusao do Na - Subsarcolema para Cleft',yaxis=dict(title= 'mM/ms',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=I_cabk_junc,type='scatter',)],layout=dict(title='Ca Background Current in Cleft',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=I_cabk_sl,type='scatter')],layout=dict(title='Ca Background Current in Subsarcollema',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=I_cabk,type='scatter')],layout=dict(title='Ca Background Current',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=I_pca_junc,type='scatter')],layout=dict(title='Ca Pump Current in Cleft',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=I_pca_sl,type='scatter',)],layout=dict(title='Ca Pump Current in Subsarcollema',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=I_pca,type='scatter',)],layout=dict(title='Ca Pump Current',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=z40,type='scatter',)],layout=dict(title='State C2 - Ca channel type L',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=z41,type='scatter',)],layout=dict(title='State C1',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=z42,type='scatter',)],layout=dict(title='State I1Ca',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=z43,type='scatter')],layout=dict(title='State I2Ca',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=z44,type='scatter')],layout=dict(title='State I1Ba',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=z45,type='scatter')],layout=dict(title='I2Ba',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=Po,type='scatter',)],layout=dict(title='State Po',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=z4041,type='scatter',name='C'),dict(x=t,y=z4243,type='scatter',name='ICa'),dict(x=t,y=z4445,type='scatter',name='IVa'),dict(x=t,y=Po,type='scatter',name='O')],layout=dict(title='LTCC',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2,),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=z0,type='scatter',)],layout=dict(title='Activation Gate - Na Fast',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),        
        dict(data=[dict(x=t,y=z1,type='scatter',)],layout=dict(title='Fast Inactivation Gate - Na Fast',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),        		
        dict(data=[dict(x=t,y=z2,type='scatter',)],layout=dict(title='Slow Inactivation Gate - Na Fast',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),                
        dict(data=[dict(x=t,y=I_tos,type='scatter',name='Itos'),dict(x=t,y=I_tof,type='scatter',name='Itof'),],layout=dict(title='Ito - K Current',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2,),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=C0fC1fC2fC3f,type='scatter',name='Fechado'),dict(x=t,y=CI0fCI1fCI2fCI3f,type='scatter',name='Fechado Inativo'),dict(x=t,y=Of,type='scatter',name='Aberto'),dict(x=t,y=If,type='scatter',name='Inativo')],layout=dict(title='Estados do Canal Itof',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2,),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=C0sC1sC2sC3s,type='scatter',name='Fechado'),dict(x=t,y=CI0sCI1sCI2sCI3s,type='scatter',name='Fechado Inativo'),dict(x=t,y=Os,type='scatter',name='Aberto'),dict(x=t,y=Is,type='scatter',name='Inativo')],layout=dict(title='Estados do Canal Itos',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2,),xaxis=dict(title= 'Time(ms)'))),		
		dict(data=[dict(x=t,y=Fcontr,type='scatter',name='Force in mN/mm^2'),dict(x=t,y=Cacy,type='scatter',name='[Ca2+]cy in mM')],layout=dict(title='Currents IKr and IKs',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),		
        dict(data=[dict(x=t,y=N0,type='scatter',)],layout=dict(title='N0',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),		
        dict(data=[dict(x=t,y=N1,type='scatter',name='N1'),dict(x=t,y=P0,type='scatter',name='P0'),dict(x=t,y=P1,type='scatter',name='P1'),dict(x=t,y=P2,type='scatter',name='P2'),dict(x=t,y=P3,type='scatter',name='P3')],layout=dict(title='Estados do Canal Itos',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2,),xaxis=dict(title= 'Time(ms)'))),		
        dict(data=[dict(x=t,y=Fcontr,type='scatter',)],layout=dict(title='Forca de Contracao',yaxis=dict(title= 'Force in mN/mm2',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=SL,type='scatter',)],layout=dict(title='Sarcomere length in um',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
        dict(data=[dict(x=t,y=CelL,type='scatter',)],layout=dict(title='Cell length in um',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),		
        dict(data=[dict(x=z37,y=Fcontr,type='scatter',)],layout=dict(title='Forca de Contracao vs Calcio Citosol',yaxis=dict(title= 'Force in mN/mm2',ticklen= 5,gridwidth= 2),xaxis=dict(title= '[Ca]cyt'))),				
        dict(data=[dict(x=t,y=IKur,type='scatter',name='IKur'),dict(x=t,y=IKss,type='scatter',name='IKss'),],layout=dict(title='Current IKur and IKss',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2,),xaxis=dict(title= 'Time(ms)')))		
        ]

    # Add "ids" to each of the graphs to pass up to the client
    # for templating
    ids = ['graph-{}'.format(i) for i, _ in enumerate(graphs)]

    # Convert the figures to JSON
    # PlotlyJSONEncoder appropriately converts pandas, datetime, etc
    # objects to their JSON equivalents
    graphJSON = json.dumps(graphs, cls=plotly.utils.PlotlyJSONEncoder)

    return render_template('resultados.html',
                           ids=ids,
                           graphJSON=graphJSON)		


@app.route('/')
def principal():
    return render_template('index.html')

@app.route('/correnteical')
def correnteical():
    return render_template('correnteical.html')

@app.route('/correntepotassio')
def correntepotassio():
    return render_template('correntepotassio.html')

@app.route('/correntesodio')
def correntesodio():
    return render_template('correntesodio.html')

@app.route('/dinamicacalcio')
def dinamicacalcio():
    return render_template('dinamicacalcio.html')

@app.route('/forcacontracao')
def forcacontracao():
    return render_template('forcacontracao.html')

@app.route('/outras')
def outras():
    return render_template('outras.html')

if __name__ == "__main__":
    app.run(debug=True)