from flask import Flask, render_template, request, url_for, flash, redirect, session, jsonify
import numpy as np
from cardiojunction import cardiojunction
import warnings
warnings.filterwarnings('ignore')
app = Flask(__name__)
app.secret_key = 'fsdfsdfsdfsdfsfsd'
import json
import plotly
import os.path
from os import path

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

    #L = L - t_ap;
	
    #L = L - t_ap;	

    Ap = int(request.args.get('txtEstimuloVoltagem'));
        
    v_resting = int(request.args.get('txtPotencialRepouso'));

    tap = int(request.args.get('txtTempoSimulacaoAplicacao'));

    w = int(request.args.get('txtLarguraPulso'));

    f = int(request.args.get('txtFrequencia'));
        
    Delay = 1000*1/f - w;

    A_inj = float(request.args.get('txtAmplitudeInjetada'));
        
    tig = int(request.args.get('txtInstantePlotagem'));
	
    arquivoNome = str(protocol) + str(Model_ICaL) + str(Model_Na) + str(Model_Ito) + str(Model_IKr) + str(Model_IKs) + str(Model_Force) + str(cellLength) + str(Lsarc) + str(BLOCKSRPUMP) + str(STIMULUSRPUMP) + str(BLOCKNCX) + str(STIMULUSNCX) + str(CAFEINA) + str(BLOCKCICR) + str(BLOCKIKs) + str(BLOCKIKr) + str(BLOCKItof) + str(BLOCKItos) + str(BLOCKINa) + str(BLOCKICaL) + str(STIMULUSICaL) + str(t_ap) + str(L) + str(Ap) + str(v_resting) + str(tap) + str(w) + str(f) + str(Delay) + str(A_inj) + str(tig);
		
    if(not path.exists(str(arquivoNome) + str(".txt"))):		

        print("");	
        c = cardiojunction();
	
        [t,I_app,v,I_Ca,I_Na,z37,z35,z36,z30,I_kr,I_ks,C0KsC1Ks,O1KsO2Ks,OKr,C1KrC2KrC3Kr,I_kp,I_ki,
        C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,O1,O2,J_SRCarel,J_serca,J_SRleak,z13,z15,z14,
        canalSerca,I_ncx,I_nak,z34,z32,z31,z33,DifNa_sl_cl,I_cabk_junc,I_cabk_sl,I_cabk,I_pca_junc,I_pca_sl,I_pca,
        Po,z40,z41,z42,z43,z44,z45,z4041,z4243,z4445,z3,z4,ONa,z484950,z515253,z5455,
        z0,z1,z2,I_tos,I_tof,Of,If,Os,Is,C0fC1fC2fC3f,CI0fCI1fCI2fCI3f,C0sC1sC2sC3s,CI0sCI1sCI2sCI3s,z7,
        z8,z9,z10,FORCA,Lsim,Cacy,CelL,Fcontr,N0,N1,P0,P1,P2,P3,SL,IKur,IKss]=c.principal(protocol, Model_ICaL, Model_Na, Model_Ito, Model_IKr, Model_IKs, Model_Force,cellLength, Lsarc, BLOCKSRPUMP,STIMULUSRPUMP,BLOCKNCX,STIMULUSNCX,CAFEINA,BLOCKCICR,BLOCKIKs,BLOCKIKr,BLOCKItof,BLOCKItos,BLOCKINa,BLOCKICaL,STIMULUSICaL,t_ap,L,Ap,v_resting,tap,w,f,Delay,A_inj,tig);

        	 
    #UM = np.linspace(1,1,len(t));
    #UM = UM.transpose();

        graphs = [
            dict(data=[dict(x=t,y=I_app,type='scatter')],layout=dict(title='Corrente Injetada',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2,),xaxis=dict(title= 'Time(ms)'))),
            dict(data=[dict(x=t,y=Fcontr,type='scatter',)],layout=dict(title='Força de Contração',yaxis=dict(title= 'Force in mN/mm2',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
            dict(data=[dict(x=t,y=SL,type='scatter',)],layout=dict(title='Sarcomere length in um',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),			
            dict(data=[dict(x=t,y=z36,type='scatter'),],layout=dict(title='CA cyt',yaxis=dict(title= '[CA2+]cyt in mM',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
            dict(data=[dict(x=t,y=SL,type='scatter',)],layout=dict(title='Sarcomere length in um',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),						
            dict(data=[dict(x=t,y=z37,type='scatter')],layout=dict(title='CA cleft',yaxis=dict(title= '[CA2+]cleft in mM',ticklen= 5,gridwidth= 2,),xaxis=dict(title= 'Time(ms)'))),			
            dict(data=[dict(x=t,y=z35,type='scatter'),],layout=dict(title='CA ss',yaxis=dict(title= '[CA2+]ss in mM',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),			
		    dict(data=[dict(x=t,y=z30,type='scatter'),],layout=dict(title='CA SR',yaxis=dict(title= '[CA2+]SR in mM',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),			
            dict(data=[dict(x=t,y=I_Ca,type='scatter'),],layout=dict(title='ICal Current',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),			
            dict(data=[dict(x=t,y=I_Na,type='scatter'),],layout=dict(title='Fast Na Current',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),			
            dict(data=[dict(x=t,y=I_ncx,type='scatter')],layout=dict(title='NCX',yaxis=dict(title= 'mM/ms',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
            dict(data=[dict(x=t,y=I_nak,type='scatter')],layout=dict(title='NaK',yaxis=dict(title= 'mM/ms',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
		    dict(data=[dict(x=t,y=I_kr,type='scatter',name='Ikr'),dict(x=t,y=I_ks,type='scatter',name='Iks')],layout=dict(title='Currents IKr and IKs',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
            dict(data=[dict(x=t,y=C0KsC1Ks,type='scatter',name='CLOSED'),dict(x=t,y=O1KsO2Ks,type='scatter',name='OPEN')],layout=dict(title='States Channel IKs',yaxis=dict(title= 'Probability',ticklen= 5,gridwidth= 2,),xaxis=dict(title= 'Time(ms)'))),
            dict(data=[dict(x=t,y=C1KrC2KrC3Kr,type='scatter',name='CLOSED'),dict(x=t,y=I_kr,type='scatter',name='INACTIVED'),dict(x=t,y=OKr,type='scatter',name='OPEN')],layout=dict(title='[CA]ss',yaxis=dict(title= 'States channel IKr',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
            dict(data=[dict(x=t,y=I_kp,type='scatter'),],layout=dict(title='IKp Current',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
            dict(data=[dict(x=t,y=I_ki,type='scatter'),],layout=dict(title='IKi Current',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
            dict(data=[dict(x=t,y=J_SRCarel,type='scatter')],layout=dict(title='CICR',yaxis=dict(title= 'mM/ms',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
            dict(data=[dict(x=t,y=J_SRleak,type='scatter')],layout=dict(title='Leak SR',yaxis=dict(title= 'mM/ms',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),            
            dict(data=[dict(x=t,y=J_serca,type='scatter')],layout=dict(title='Serca Pump SR',yaxis=dict(title= 'mM/ms',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)')))		
            ]
 

    # Add "ids" to each of the graphs to pass up to the client
    # for templating
        ids = ['graph-{}'.format(i) for i, _ in enumerate(graphs)]

    # Convert the figures to JSON
    # PlotlyJSONEncoder appropriately converts pandas, datetime, etc
    # objects to their JSON equivalents
        graphJSON = json.dumps(graphs, cls=plotly.utils.PlotlyJSONEncoder);
        arquivo = open(str(arquivoNome) + str(".txt"),'w');
        arquivo.write(graphJSON);
        arquivo.close();
        print("");				
    else:	
        print("");
        arquivo = open(str(arquivoNome) + str(".txt"),'r');
        graphJSON =arquivo.read();        
        obj = json.loads(graphJSON, cls=json.JSONDecoder)		
        ids = ['graph-{}'.format(i) for i, _ in enumerate(obj)]			
        arquivo.close();
        print("");		
		
    return render_template('resultados.html',
                           ids=ids,
                           graphJSON=graphJSON)		
    #return render_template('resultado.html')

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
    app.run(debug=True, host='0.0.0.0')
#    app.run(host='0.0.0.0')
#    app.run(debug=True)
