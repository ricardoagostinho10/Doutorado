$( document ).ready(function() {
		pega_dados();
	});

	window.setTimeout(function() {
		$(".alert").fadeTo(500, 0).slideUp(500, function(){
			$(this).remove();
		});
	}, 3000)

	var selecionado = {
		nome : "",
		dado : [],
		grandeza : "",
		subtitutlo : ""
	};
	var inicio = true;

	var t;
	var I_app;


	function gera_grafico(){
		Highcharts.chart('container', {
		chart: {
				zoomType: 'x'
			},
		title: {
			text: selecionado.nome
		},
		subtitle: {
			text: selecionado.titulo
		},
		yAxis: {
			title: {
				text: selecionado.grandeza
			}
		},
		legend: {
			layout: 'vertical',
			align: 'right',
			verticalAlign: 'middle'
		},
		plotOptions: {
			series: {
				pointStart: 0,
				animation: {
					duration: 3000
					//easing: 'easeOutBounce'
				}
			}
		},
		series: [{
			name: selecionado.grandeza,
			data: selecionado.dado,
			type: 'spline',
			}]
		})
	}

	function pega_dados(){
		selProtocolo = $('#selProtocolo').val();
			
		selCalcio = $('#selCalcio').val();

		selSodio = $('#selSodio').val();

		selPotassioTrans =  $('#selPotassioTrans').val();

		selPotassioR = $('#selPotassioR').val();

		selPotassioS =  $('#selPotassioS').val();

		selForca = $('#selForca').val();

		txtComCel  = $('#txtComCel').val();

		txtComSar = $('#txtComSar').val();
		
		txtBloqueioPump = $('#txtBloqueioPump').val();
		txtEstimuloPump = $('#txtEstimuloPump').val();
		txtBloqueioNCX = $('#txtBloqueioNCX').val();
		txtEstimuloNCX = $('#txtEstimuloNCX').val();
		selCafeina = $('#selCafeina').val();
		txtBloqueioCICR = $('#txtBloqueioCICR').val();
		txtBloqueioIKs = $('#txtBloqueioIKs').val();
		txtBloqueioIKr = $('#txtBloqueioIKr').val();
		txtBloqueioIKtof = $('#txtBloqueioIKtof').val();
		txtBloqueioIKtos = $('#txtBloqueioIKtos').val();
		txtBloqueioINaFast = $('#txtBloqueioINaFast').val();
		txtBloqueioICaL = $('#txtBloqueioICaL').val();
		txtEstimuladorICaL = $('#txtEstimuladorICaL').val();
	   
		txtInstanteAplicacaoPulso = $('#txtInstanteAplicacaoPulso').val();
		   
		txtTempoSimulacao = $('#txtTempoSimulacao').val();


		txtEstimuloVoltagem = $('#txtEstimuloVoltagem').val();
			
		txtPotencialRepouso = $('#txtPotencialRepouso').val();

		txtTempoSimulacaoAplicacao = $('#txtTempoSimulacaoAplicacao').val();

		txtLarguraPulso = $('#txtLarguraPulso').val();

		txtFrequencia = $('#txtFrequencia').val();
			
		txtAmplitudeInjetada = $('#txtAmplitudeInjetada').val();
			
		txtInstantePlotagem = $('#txtInstantePlotagem').val();

			
		if (!selProtocolo)
		{
			selProtocolo = 1;
		}
		if (!selCalcio)
		{
			selCalcio = 1;
		}
		if (!selSodio)
		{
			selSodio = 0;
		}
		if (!selPotassioTrans)
		{
			selPotassioTrans = 1;
		}
		if (!selPotassioR)
		{
			selPotassioR = 1;
		}
		if (!selPotassioS)
		{
			selPotassioS = 1;
		}
		if (!selForca)
		{
			selForca = 1;
		}	
		if (!txtComCel)
		{
			txtComCel = 100;
		}		
		if (!txtComSar)
		{
			txtComSar = 1.05;
		}	
		if (!txtBloqueioPump)
		{
			txtBloqueioPump = 0;
		}	
		if (!txtEstimuloPump)
		{
			txtEstimuloPump = 0;
		}	
		if (!txtBloqueioNCX)
		{
			txtBloqueioNCX = 0;
		}    
		if (!txtEstimuloNCX)
		{
			txtEstimuloNCX = 0;
		}  
		if (!selCafeina)
		{
			selCafeina = 0;
		} 
		if (!txtBloqueioCICR)
		{
			txtBloqueioCICR = 0;
		} 
		if (!txtBloqueioIKs)
		{
			txtBloqueioIKs = 0;
		} 
		if (!txtBloqueioIKr)
		{
			txtBloqueioIKr = 0;
		} 
		if (!txtBloqueioIKtof)
		{
			txtBloqueioIKtof = 0;
		} 
		if (!txtBloqueioIKtos)
		{
			txtBloqueioIKtos = 1.05;
		} 
		if (!txtBloqueioINaFast)
		{
			txtBloqueioINaFast = 0;
		} 
		if (!txtBloqueioICaL)
		{
			txtBloqueioICaL = 0;
		} 
		if (!txtEstimuladorICaL)
		{
			txtEstimuladorICaL = 0;
		} 	
		if (!txtInstanteAplicacaoPulso)
		{
			txtInstanteAplicacaoPulso = 1000;
		} 
		if (!txtTempoSimulacao)
		{
			txtTempoSimulacao = 5000;
		} 
		if (!txtEstimuloVoltagem)
		{
			txtEstimuloVoltagem = 15;
		} 
		if (!txtPotencialRepouso)
		{
			txtPotencialRepouso = -80;
		} 
		if (!txtTempoSimulacaoAplicacao)
		{
			txtTempoSimulacaoAplicacao = 0;
		}
		if (!txtLarguraPulso)
		{
			txtLarguraPulso = 200;
		} 
		if (!txtFrequencia)
		{
			txtFrequencia = 1;
		} 
		if (!txtAmplitudeInjetada)
		{
			txtAmplitudeInjetada = 9.5;
		} 
		if (!txtInstantePlotagem)
		{
			txtInstantePlotagem = 0;
		}
		
		$.getJSON('/dado_grafico.json',
		{
		selProtocolo:selProtocolo,        
		selCalcio:selCalcio,
		selSodio:selSodio,
		selPotassioTrans:selPotassioTrans,
		selPotassioR:selPotassioR,
		selPotassioS:selPotassioS,
		selForca:selForca,
		txtComCel:txtComCel,
		txtComSar: txtComSar,
		txtBloqueioPump:txtBloqueioPump,
		txtEstimuloPump:txtEstimuloPump ,
		txtBloqueioNCX:txtBloqueioNCX ,
		txtEstimuloNCX:txtEstimuloNCX ,
		selCafeina:selCafeina ,
		txtBloqueioCICR:txtBloqueioCICR ,
		txtBloqueioIKs:txtBloqueioIKs ,
		txtBloqueioIKr:txtBloqueioIKr ,
		txtBloqueioIKtof:txtBloqueioIKtof ,
		txtBloqueioIKtos:txtBloqueioIKtos ,
		txtBloqueioINaFast:txtBloqueioINaFast ,
		txtBloqueioICaL:txtBloqueioICaL ,
		txtEstimuladorICaL:txtEstimuladorICaL ,   
		txtInstanteAplicacaoPulso:txtInstanteAplicacaoPulso ,       
		txtTempoSimulacao:txtTempoSimulacao ,
		txtEstimuloVoltagem:txtEstimuloVoltagem ,        
		txtPotencialRepouso:txtPotencialRepouso ,
		txtTempoSimulacaoAplicacao:txtTempoSimulacaoAplicacao ,
		txtLarguraPulso:txtLarguraPulso ,
		txtFrequencia:txtFrequencia ,        
		txtAmplitudeInjetada:txtAmplitudeInjetada ,        
		txtInstantePlotagem:txtInstantePlotagem	
		})
		.done(function( data ){
			console.log("resposta");
			t = data[0];
			I_app = data[1];
			if (inicio)
			{
				selecionado.dado = I_app;
				selecionado.nome = "I_app";
				selecionado.grandeza = "Hz";
				selecionado.titulo = "M4 - sei la o que é isso";
			}
			inicio = false;
			gera_grafico();
		});
	}

	$(".dropdown-menu li a").click(function(){
		$(this).parents(".dropdown").find('.btn').html($(this).text() + ' <span class="caret"></span>');
		$(this).parents(".dropdown").find('.btn').val($(this).data('value'));
		var sel = $(this).text();
		switch (sel) {
			case 'M4':
				selecionado.dado = I_app;
				selecionado.nome = "M4";
				selecionado.grandeza = "Hz";
				selecionado.titulo = "M4 - sei la o que é isso";
				break;
		default:
				selecionado.dado = I_app;
				selecionado.nome = "M4";
				selecionado.grandeza = "Hz";
				selecionado.titulo = "M4 - sei la o que é isso";
	  }
	  gera_grafico();
	});

	$("#Atualizar").click(function(){
		gera_grafico();
	});

	$("#recalcular").click(function(){
		pega_dados();
	});