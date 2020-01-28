from flask import Flask, render_template, request, url_for, flash, redirect, session, jsonify
import numpy as np
from cardiojunction import cardiojunction
import warnings
warnings.filterwarnings('ignore')
app = Flask(__name__)
app.secret_key = 'fsdfsdfsdfsdfsfsd'
import json
import plotly

#import os.path
#from os import path

@app.route('/graficos', methods = ['POST'])
def get_plot():    

	
    protocol = int(request.form['selProtocolo']);
        
    Model_ICaL = int(request.form['selCalcio']);

    Model_Na = int(request.form['selSodio']);

    Model_Ito =  int(request.form['selPotassioTrans']);

    Model_IKr = int(request.form['selPotassioR']);

    Model_IKs =  int(request.form['selPotassioS']);

    Model_Force = int(request.form['selForca']);

    cellLength  = int(request.form['txtComCel']);

    Lsarc = float(request.form['txtComSar']);
	
    BLOCKSRPUMP = float(request.form['txtBloqueioPump'])/100;
    STIMULUSRPUMP = float(request.form['txtEstimuloPump'])/100;
    BLOCKNCX = float(request.form['txtBloqueioNCX'])/100;
    STIMULUSNCX = float(request.form['txtEstimuloNCX'])/100;
    CAFEINA = int(request.form['selCafeina'])/100;
    BLOCKCICR = float(request.form['txtBloqueioCICR'])/100;
    BLOCKIKs = float(request.form['txtBloqueioIKs'])/100;
    BLOCKIKr = float(request.form['txtBloqueioIKr'])/100;
    BLOCKItof = float(request.form['txtBloqueioIKtof'])/100;
    BLOCKItos = float(request.form['txtBloqueioIKtos'])/100;
    BLOCKINa = float(request.form['txtBloqueioINaFast'])/100;
    BLOCKICaL = float(request.form['txtBloqueioICaL'])/100;
    STIMULUSICaL = float(request.form['txtEstimuladorICaL'])/100;

    if(protocol==1):
        	
        t_ap = int(request.form['txtInstanteAplicacaoPulsoPotencial']);
       
        L = int(request.form['txtTempoSimulacaoPotencial']);

        Ap = int(request.form['txtEstimuloVoltagemPotencial']);
        
        v_resting = int(request.form['txtPotencialRepousoPotencial']);

        tap = int(request.form['txtTempoSimulacaoAplicacaoPotencial']);

        w = int(request.form['txtLarguraPulsoPotencial']);
#        if(w<10):
#            w = 10;		
		
        f = int(request.form['txtFrequenciaPotencial']);
        
        Delay = 1000*1/f - w;

        A_inj = float(request.form['txtAmplitudeInjetadaInjetada']);
        
        tig = int(request.form['txtInstantePlotagemPotencial']);
    else:
        t_ap = int(request.form['txtInstanteAplicacaoPulsoInjetada']);
       
        L = int(request.form['txtTempoSimulacaoInjetada']);

        Ap = int(request.form['txtEstimuloVoltagemPotencial']);
        
        v_resting = int(request.form['txtPotencialRepousoPotencial']);

        tap = int(request.form['txtTempoSimulacaoAplicacaoInjetada']);

        w = int(request.form['txtLarguraPulsoInjetada']);
#        if(w<10):
#            w = 10;		
		
        f = int(request.form['txtFrequenciaInjetada']);
        
        Delay = 1000*1/f - w;

        tig = int(request.form['txtInstantePlotagemInjetada']);
		
        A_inj = float(request.form['txtAmplitudeInjetadaInjetada']);

    print(protocol);	    
    print(t_ap);
    print(L);
    print(Ap);
    print(v_resting);
    print(tap);
    print(f);    
    print(Delay);
    print(tig);
    print(A_inj);	
    print(w);
    c = cardiojunction();
	
    [t,I_app,v,I_Ca,I_Na,z37,z35,z36,z30,I_kr,I_ks,C0KsC1Ks,O1KsO2Ks,OKr,C1KrC2KrC3Kr,I_kp,I_ki,
    C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,O1,O2,J_SRCarel,J_serca,J_SRleak,z13,z15,z14,
    canalSerca,I_ncx,I_nak,z34,z32,z31,z33,DifNa_sl_cl,I_cabk_junc,I_cabk_sl,I_cabk,I_pca_junc,I_pca_sl,I_pca,
    Po,z40,z41,z42,z43,z44,z45,z4041,z4243,z4445,z3,z4,ONa,z484950,z515253,z5455,
    z0,z1,z2,I_tos,I_tof,Of,If,Os,Is,C0fC1fC2fC3f,CI0fCI1fCI2fCI3f,C0sC1sC2sC3s,CI0sCI1sCI2sCI3s,z7,
    z8,z9,z10,FORCA,Lsim,Cacy,CelL,Fcontr,N0,N1,P0,P1,P2,P3,SL,IKur,IKss]=c.principal(protocol, Model_ICaL, Model_Na, Model_Ito, Model_IKr, Model_IKs, Model_Force,cellLength, Lsarc, BLOCKSRPUMP,STIMULUSRPUMP,BLOCKNCX,STIMULUSNCX,CAFEINA,BLOCKCICR,BLOCKIKs,BLOCKIKr,BLOCKItof,BLOCKItos,BLOCKINa,BLOCKICaL,STIMULUSICaL,t_ap,L,Ap,v_resting,tap,w,f,Delay,A_inj,tig);

    graphs = [
         dict(data=[dict(x=t,y=v,type='scatter')],layout=dict(title='Voltagem',yaxis=dict(title= 'm/V',ticklen= 5,gridwidth= 2,),xaxis=dict(title= 'Time(ms)'))),
         dict(data=[dict(x=t,y=z36,type='scatter'),],layout=dict(title='Ca<sup>2+</sup>CYT',yaxis=dict(title= '[Ca<sup>2+</sup>]cyt in mM',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
         dict(data=[dict(x=t,y=Fcontr,type='scatter',)],layout=dict(title='Força de Contração',yaxis=dict(title= 'Force in mN/mm<sup>2</sup>',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
         dict(data=[dict(x=t,y=SL,type='scatter',)],layout=dict(title='SL',yaxis=dict(title= 'Sarcomere length in um',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),			
         dict(data=[dict(x=t,y=z35,type='scatter')],layout=dict(title='Ca<sup>2+</sup>Cleft',yaxis=dict(title= '[Ca<sup>2+</sup>]cleft in mM',ticklen= 5,gridwidth= 2,),xaxis=dict(title= 'Time(ms)'))),			
         dict(data=[dict(x=t,y=z36,type='scatter'),],layout=dict(title='Ca<sup>2+</sup>SS',yaxis=dict(title= '[Ca<sup>2+</sup>]ss in mM',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),			
         dict(data=[dict(x=t,y=z37,type='scatter'),],layout=dict(title='Ca<sup>2+</sup>CY',yaxis=dict(title= '[Ca<sup>2+</sup>]cyt in mM',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),		 
	     dict(data=[dict(x=t,y=z30,type='scatter'),],layout=dict(title='Ca<sup>2+</sup>SR',yaxis=dict(title= '[Ca<sup>2+</sup>]SR in mM',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),			
         dict(data=[dict(x=t,y=I_Ca,type='scatter'),],layout=dict(title='ICal Current',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),			
         dict(data=[dict(x=t,y=I_Na,type='scatter'),],layout=dict(title='INa Current',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),			
         dict(data=[dict(x=t,y=I_ncx,type='scatter')],layout=dict(title='INaCa',yaxis=dict(title= 'mM/ms',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
         dict(data=[dict(x=t,y=I_nak,type='scatter')],layout=dict(title='NaK',yaxis=dict(title= 'mM/ms',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
         dict(data=[dict(x=t,y=J_SRCarel,type='scatter')],layout=dict(title='CICR',yaxis=dict(title= 'mM/ms',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
         dict(data=[dict(x=t,y=J_SRleak,type='scatter')],layout=dict(title='Leak SR',yaxis=dict(title= 'mM/ms',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),            
         dict(data=[dict(x=t,y=J_serca,type='scatter')],layout=dict(title='Serca Pump SR',yaxis=dict(title= 'mM/ms',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),
         dict(data=[dict(x=t,y=I_pca_sl,type='scatter')],layout=dict(title='Serca Pump SL',yaxis=dict(title= 'mM/ms',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),		 
         dict(data=[dict(x=t,y=I_tof,type='scatter')],layout=dict(title='IKtof',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),		 
         dict(data=[dict(x=t,y=I_ki,type='scatter'),],layout=dict(title='IK<sub>i</sub> Current',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),		 
         dict(data=[dict(x=t,y=I_ks,type='scatter')],layout=dict(title='Iks',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),		 		 
         dict(data=[dict(x=t,y=I_kr,type='scatter')],layout=dict(title='Ikr',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),		 		 		 
         dict(data=[dict(x=t,y=IKur,type='scatter')],layout=dict(title='IKur',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)'))),		 		 		          
         dict(data=[dict(x=t,y=IKss,type='scatter')],layout=dict(title='IKss',yaxis=dict(title= 'uA/uF',ticklen= 5,gridwidth= 2),xaxis=dict(title= 'Time(ms)')))		 		 		          		 
         ]
 

    # Add "ids" to each of the graphs to pass up to the client
    # for templating
    ids = ['graph-{}'.format(i) for i, _ in enumerate(graphs)]

    # Convert the figures to JSON
    # PlotlyJSONEncoder appropriately converts pandas, datetime, etc
    # objects to their JSON equivalents
    graphJSON = json.dumps(graphs, cls=plotly.utils.PlotlyJSONEncoder);
    	
    return render_template('resultados.html',
                           ids=ids,
                           graphJSON=graphJSON)		
    #return render_template('resultado.html')

@app.route('/')
def principal():
    return render_template('index.html')

@app.route('/parametros')
def parametros():
    return render_template('parametros.html')

@app.route('/resultados')
def resultados():
    return render_template('resultados.html')

@app.route('/sobre')
def sobre():
    return render_template('sobre.html')
	
if __name__ == "__main__":
#    try:
#        import psyco
#        print("opa");		
#        psyco.full();
#    except: print('Psyco optimizer not installed');
    
    #app.run(debug=True)
    app.run(debug=False, host='0.0.0.0')
    
    #app.run(debug=False, host='0.0.0.0')
