pkg load signal;
close all;

%% 1. Cargar el archivo intensidad_RGB.mat y visualizar la señal RGB en
% superposición (3 canales). Elija alguno de los 3 canales que, según su criterio, posea
% la señal más útil a los efectos de analizar la dinámica de los latidos y utilícelo para
% realizar los puntos siguientes.
%
printf('\nEjercicio 1:\n');

% Cargo los datos.
printf('Cargando los datos.\n')
intensity = load('intensidad_RGB.mat');

% Las señales.
printf('Cargando las señales.\n')
rgb = intensity.brillo;
%rgb = intensity.intensidad(1:2200, :);
redSignal = rgb(:, 1);
greenSignal = rgb(:, 2);
blueSignal = rgb(:, 3);
N = length(rgb)

% La frecuencia.
printf('Cargando la frecuencia.\n')
samplingPeriod = 1/intensity.FR;

% Calculo el eje temporal.
printf('Calculando el eje temporal.\n')
timeAxis = [0 : samplingPeriod : samplingPeriod * (N-1)];

% Normalizo las señales para que tengan valor máximo unitario.
printf('Normalizando las señales.\n')
redSignal = redSignal / max(redSignal);
greenSignal = greenSignal / max(greenSignal);
blueSignal = blueSignal / max(blueSignal);

% Dibujo en superposición las señales.
printf('Dibujando en superposición.\n')
figure;
plot(timeAxis, redSignal, 'r');
grid on;
set(gca, 'xminorgrid', 'on');
hold on;
plot(timeAxis, greenSignal, 'g');
hold on;
plot(timeAxis, blueSignal, 'b');
xlabel('Tiempo [s]');
ylabel('Intensidad');
legend('Roja', 'Verde', 'Azul');
title('Intensidad por color.')

% Guardo en un archivo.
print -dpng 'images/ej1/superposicion.png';

% Elijo una señal para el resto de los items.
signal = greenSignal;





%% 2. Del gráfico del punto anterior, estimar aproximadamente los latidos por minuto
% (LPM). Identifique el momento a partir del cual la frecuencia cardiaca comienza a
% incrementarse, y en cuánto se incrementa.
%
printf('\nEjercicio 2:\n');

% Valor mínimo de la señal.
minimum = min(signal);

% Busco los picos de la señal con el cero en 'minimum'.
printf('Buscando los picos en la señal.\n');
[signalPeaksValues signalPeakFrames c] = findpeaks(signal .- minimum, 'MinPeakHeight', 0, 'MinPeakDistance', 10);

% Corrijo el gap en los valores de los picos.
signalPeaksValues = signalPeaksValues .+ minimum;

% Estimo el promedio de latidos por minuto.
rate = length(signalPeaksValues)*60/(max(timeAxis));
printf('Promedio de latidos por minuto estimado: %d LPM.\n', rate);

% Este grafico muestra una tendencia decreciente de
% los períodos (una frecuencia creciente) desde 1.1 hasta 0.73, 
% digamos: 1/freq = 0.915 (65.5 LPM) > 1Hz
%
% Luego, en aproximadamente el instante t = 50s,
% se estanca en, digamos: 1/freq = 0.63 (95 LPM) < 5Hz
%
printf('Mostrando los períodos.\n');
figure;
deltas = diff(timeAxis(signalPeakFrames));
deltas = [0, deltas];
stem(timeAxis(signalPeakFrames), deltas);
grid on;
set(gca, 'xminorgrid', 'on');
xlabel('Tiempo [s]');
ylabel('Intervalos entre latidos [s]');
title('Evolución de los intervalos entre latidos.');
print -dpng 'images/ej2/intervalos.png';

% Las frecuencias estimadas usando los datos anteriores
printf('Mostrando las frecuencias instantáneas.\n');
figure;
stem(timeAxis(signalPeakFrames), 60 ./ deltas);
grid on;
set(gca, 'xminorgrid', 'on');
xlabel('Tiempo [s]');
ylabel('Frecuencia cardíaca [LPM]');
title('Evolución de las frecuencias instantáneas estimadas.');
print -dpng 'images/ej2/frecuencias.png';

% Los picos en superposición con la señal.
figure;
printf('Mostrando los picos en superposición con la señal.\n');
plot(timeAxis, signal, 'r');
grid on, hold on;
set(gca, 'xminorgrid', 'on');
stem(timeAxis(signalPeakFrames), signalPeaksValues);
xlabel('Tiempo [s]');
ylabel('Intensidad');
legend('Señal', 'Picos detectados');
title('Picos en la señal.')
print -dpng 'images/ej2/picos.png';





%% 3. Rehacer el punto anterior, pero utilizando DFT. ¿Puede en este caso identificar
% cuándo se da el cambio de LPM?
%
printf('\nEjercicio 3:\n');

% Calculo el eje de frecuencia discreta.
discreteFrequencyStep = 2*pi / N;
discreteFrequencyAxis = [-pi : discreteFrequencyStep : pi-discreteFrequencyStep];

% Calculo el eje de frecuencia angular y en Hertz.
continuousFrequencyAxis = discreteFrequencyAxis / samplingPeriod;
hertzAxis = continuousFrequencyAxis / (2*pi);

% Calculo el espectro y lo ubico entre -pi y pi.
printf('Calculando la DFT.\n');
dft = fftshift(fft(signal));

function plotDFTWithXlims (hertzAxis, dft, xlims)
	figure;

	subplot(2, 1, 1);
	absolute = abs(dft);
	plot(hertzAxis, absolute);
	set(gca, 'xminorgrid', 'on');
	xlabel('Frecuencia [Hz]');
	ylabel('Modulo de la DFT');

	if length(xlims) > 0
		xlim(xlims);
	end

	title('Transformada discreta de Fourier');

	subplot(2, 1, 2);
	argument = arg(dft);
	plot(hertzAxis, argument);
	xlabel('Frecuencia [Hz]');
	ylabel('Fase de la DFT');

endfunction

function plotDFT (hertzAxis, dft)
	plotDFTWithXlims(hertzAxis, dft, []);
endfunction


% En este grafico se ve un pico en unos 1.6Hz (96 LPM),
% También se ven varios picos parecidos entre 40 y 80 LPM,
% que se condicen con la frecuencia creciente que veíamos antes de los
% 50 segundos.
%
% El problema de este gráfico es que no nos da información sobre cuándo
% se ven esas frecuencias (no tiene localización temporal).
%
plotDFTWithXlims(hertzAxis, dft, [-2 2]);
print -dpng 'images/ej3/dft.png';





%% 4. Diseñar un filtro pasa-banda tipo Butterworth con banda de paso entre 0.5 Hz y 10
% Hz. Graficar respuesta en frecuencia (módulo y fase), diagrama de polos y ceros y
% respuesta al impulso (sugerencia: usar las funciones butter y fvtool de Matlab).
% ¿Qué papel juega el orden del filtro seleccionado en su diseño?
%
printf('\nEjercicio 4:\n');

function [bButter aButter] = Butterworth (filterOrder, samplingPeriod, N, timeAxis, mustPlot)
	printf('Usando Butterworth de orden %d:\n', filterOrder);

	% Paso a frecuencia discreta
	lowDiscreteFrequency = 0.5 * 2*pi * samplingPeriod;
	highDiscreteFrequency = 10 * 2*pi * samplingPeriod;

	% Armo el filtro
	printf('Calculando los coeficientes para el Butterworth.\n');
	[bButter, aButter] = butter(filterOrder/2, [lowDiscreteFrequency/pi, highDiscreteFrequency/pi]);

	% Calculo la respuesta en frecuencia
	printf('Calculando la respuesta en frecuencia del Butterworth.\n');
	[HButter, OmegaButter] = freqz(bButter, aButter);

	% Tomo las señales y ejes a graficar
	absolute = abs(HButter);
	phase = arg(HButter);
	OmegaButter = OmegaButter/pi;

	% Grafico la respuesta en frecuencia
	if mustPlot
		figure;

		subplot(2, 1, 1);
		plot(OmegaButter, absolute);
		set(gca, 'xminorgrid', 'on');
		xlabel('Frecuencia normalizada');
		ylabel('Modulo');
		ylim([0 1]);
		title('Respuesta en frecuencia');

		subplot(2, 1, 2);
		plot(OmegaButter, phase);
		set(gca, 'xminorgrid', 'on');
		% Limito a la parte que me interesa
		%xlim([lowDiscreteFrequency/pi highDiscreteFrequency/pi]);
		xlabel('Frecuencia normalizada');
		ylabel('Fase');

		fileName = strcat('images/ej4/respuestafrecuencia', int2str(filterOrder), '.png');
		print(fileName, '-dpng');

		% Diagrama de polos y ceros
		figure;
		zplane(bButter, aButter);
		title('Diagrama de polos y ceros');
		legend('Circunferencia unitaria', 'Polo', 'Cero');

		fileName = strcat('images/ej4/polosyceros', int2str(filterOrder), '.png');
		print(fileName, '-dpng');

		% Calculo la respuesta al impulso
		printf('Calculando la respuesta al impulso del Butterworth.\n');
		impulseResponse = impz(bButter, aButter, N, 1/samplingPeriod);

		% Grafico la respuesta al impulso
		figure;
		plot(timeAxis, impulseResponse);
		set(gca, 'xminorgrid', 'on');
		xlabel('Tiempo [s]');
		ylabel('Respuesta al impulso');
		xlim([0 8]);
		ylim([-0.8 0.8]);
		title('Respuesta al impulso');
		fileName = strcat('images/ej4/respuestaimpulso', int2str(filterOrder), '.png');
		print(fileName, '-dpng');
	end

endfunction

% Vemos que las respuestas en frecuencia son más pronunciadas y la fase es sensible
% a medida que se sube el orden de los filtros. También se ve que la respuesta al
% impulso tiene un transitorio cada vez más largo.
%
Butterworth(2, samplingPeriod, N, timeAxis, true);
Butterworth(4, samplingPeriod, N, timeAxis, true);
Butterworth(6, samplingPeriod, N, timeAxis, true);
[bButter aButter] = Butterworth(10, samplingPeriod, N, timeAxis, true);





%% 5. Filtrar la señal FPG utilizando el filtro diseñado en el punto anterior mediante la
% función filter. Grafique en superposición la señal original con la filtrada y
% comente acerca de:
%	a. Remoción de derivas
%	b. Cambios en la forma de la señal
%	c. Retardo de la señal filtrada respecto de la original
%
printf('\nEjercicio 5:\n');

% Filtro la señal
printf('Filtrando la señal.\n');
butterworthFilteredSignal = filter(bButter, aButter, signal);

printf('Graficando la señal original y filtrada.\n');
figure;

% Dibujo la señal original.
plot(timeAxis, signal, 'r');
hold on;

% Dibujo la señal filtrada con el Butterworth:
%	a. Se ve que desaparece la señal de frecuencia baja (sobre la que
%	parece montada la señal original), y que desaparecen los picos
%	más bruscos, ya que se están filtrando las frecuencias altas.
%	b. Pero también se observa una deformación de la señal, dada por
%	lo mucho que varía la fase.
%	c. Se ve una diferencia de fases de 4 muestras.
%
plot(timeAxis, -butterworthFilteredSignal, 'g');
xlabel('Tiempo [s]');
ylabel('Respuesta al filtro');
grid on;
set(gca, 'xminorgrid', 'on');
title('Señal filtrada con un Butterworth');
legend('Señal original', 'Señal filtrada');
print -dpng 'images/ej5/superposicion.png';

figure;
plot(signal, 'r');
hold on;
plot(-butterworthFilteredSignal, 'g');
xlabel('Número de muestra');
ylabel('Respuesta al filtro');
grid on;
set(gca, 'xminorgrid', 'on');
title('Señal filtrada con un Butterworth');
legend('Señal original', 'Señal filtrada');
ylim([-0.05 0.05]);
xlim([1690 1710]);

print -dpng 'images/ej5/zoom.png';





%% 6. A partir de la respuesta en fase del filtro, calcule su retardo temporal 
% y compare con lo observado en el punto 5c.
%
% En la figura se ve que para las frecuencias altas, la fase es de aproximadamente
% unas cuatro muestras, lo que se condice con lo observado en el punto 5-c, ya que
% para estimar ese valor se habían comparado las posiciones de los picos
%
printf('\nEjercicio 6:\n');

[HButter, OmegaButter] = freqz(bButter, aButter);
phase = arg(HButter);
drift = (phase(3:end) - phase(1:end-2)) ./ (OmegaButter(3) - OmegaButter(1));
OmegaButter = OmegaButter(2:end-1);

figure;
plot(OmegaButter/pi, -drift);
xlabel('Frecuencia normalizada [pi x rad/s]');
ylabel('Retardo de grupo [muestras]');
grid on;
set(gca, 'xminorgrid', 'on');
title('Retardo temporal');
ylim([0 5]);

print -dpng 'images/ej6/derivada.png';





%% 7. Implementar un filtrado IIR ida y vuelta para anular la fase del filtro 
% (puede utilizar la función filtfilt de Matlab). Justificar teóricamente el 
% funcionamiento de este tipo de filtrado y cuál resulta su ventaja. Filtrar 
% nuevamente la señal FPG y comparar el resultado con lo obtenido en el punto
% anterior, particularmente en la forma de la señal y su retardo.
%
printf('\nEjercicio 7:\n');

% Señal filtrada con el IIR (filtfilt).
%	A diferencia de la señal filtrada con el Butterworth, esta señal
%	conserva la forma de la señal original, porque H x H* == |H|^2,
%	que es real (su fase es igual a cero). Que la fase sea igual a cero
%	también hace que la señal de salida no tenga retardos.
%
printf('Calculando la respuesta al filtro ida y vuelta.\n');
iirFilteredSignal = filtfilt(bButter, aButter, signal);

figure;
plot(timeAxis, iirFilteredSignal, 'b');
xlabel('Tiempo [s]');
ylabel('Respuesta');
grid on;
set(gca, 'xminorgrid', 'on');
title('Respuesta al filtro ida y vuelta');
print -dpng 'images/ej7/respuestaiir.png';





%% 8. Realizar un espectrograma de la señal antes y después de filtrar, mediante
% la función spectrogram de Matlab (sugerencia: utilice la función caxis para
% saturar los colores del espectrograma y lograr una mejor visualización).
% Justificar la longitud de ventana elegida y comente acerca del resultado
% obtenido, relacionándolo con los puntos 2 y 3. Calcule la resolución en
% frecuencia de la ventana mediante DFT en Matlab. ¿Cómo haría para obtener mejor
% resolución en frecuencia y qué se pierde con esto?
%
printf('\nEjercicio 8:\n');

function generateSpectrograms (signal, iirFilteredSignal, samplingPeriod, specgramWindow, windowWidth)
	figure;
	subplot(2, 1, 1);
	specgram(signal, windowWidth, 1/samplingPeriod, specgramWindow, round(3*windowWidth/4));
	ylim([0 4]);
	colormap('jet');
	xlim([0 80]);
	%caxis([-25 20]);
	caxis([-25 3]);
	colorbar();
	title('Señal original');

	subplot(2, 1, 2);
	specgram(iirFilteredSignal, windowWidth, 1/samplingPeriod, specgramWindow, round(3*windowWidth/4));
	ylim([0 4]);
	xlim([0 80]);
	colormap('jet');
	colorbar();
	%caxis([-25 15]);
	caxis([-25 3]);
	title('Señal filtrada');
endfunction

% Muestro DFT de la señal filtrada
dftFiltered = fftshift(fft(iirFilteredSignal));
plotDFT(hertzAxis, dftFiltered);
print -dpng 'images/ej8/dftfiltrada.png';

% Calculo el ancho de la ventana a partir de lo observado en la DFT
windowWidth = ceil((1+2/0.4)/samplingPeriod);
printf('Calculo el ancho de la ventana a partir de lo observado en la DFT: %d.\n', windowWidth);

% En el espectrograma se ve mucho ruido al principio, y luego una frecuencia creciente
% que arranca en 1Hz aproximadamente a los 20 segundos, y después de los 50 segundos
% se estanca en unos 1.6Hz.
%
% La diferencia entre la señal original y la filtrada se ve en las frecuencias bajas,
% que en la filtrada tienen un color mucho más frío.
%
% La diferencia entre la ventana de Hamming y la de Hanning se ve en el ancho/contraste
% de las curvas que se forman: la ventana de Hamming genera curvas de mayor contraste,
% mientras que las curvas generadas por la de Hanning son más gruesas.
%
specgramWindow = window(@hanning, windowWidth);
printf('Generando espectrograma con ventana Hanning\n');
generateSpectrograms(signal, iirFilteredSignal, samplingPeriod, specgramWindow, windowWidth);
print -dpng 'images/ej8/hanning.png';

printf('Generando espectrograma con ventana Hamming\n');
specgramWindow = window(@hamming, windowWidth);
generateSpectrograms(signal, iirFilteredSignal, samplingPeriod, specgramWindow, windowWidth);
print -dpng 'images/ej8/hamming.png';

% Para ganar resolución en frecuencia se puede agrandar el tamaño de la ventana,
% pero eso hace que perdamos localización temporal.





%% 9. Identificar en el espectrograma la zona donde el pulso se acelera. Observar
% con detenimiento los componentes de frecuencia que posee la señal y justificar
% el origen de cada uno (para esto último, necesitará hacer uso de la señal
% audio_det.mat para explicar todos los componentes observados).
printf('\nEjercicio 9:\n');

% Las observaciones del espectrograma están en el ejercicio anterior. La frecuencia
% que ronda 1Hz debería ser la relacionada con el pulso cardíaco, y la de 2~3Hz es
% un armónico, la otra entonces debería estar relacionada con la señal de audio...

% Cargo la señal de audio
printf('Cargando señal de audio.\n');
audioMat = load('audio_det.mat');
audioSamplingPeriod = 1/audioMat.fsa;
audio = audioMat.xa_det;
audioN = length(audio);
audioTimeAxis = [0 : audioSamplingPeriod : audioSamplingPeriod * (audioN - 1)];

% Normalizo
audio = audio / max(audio);

[bButterAudio aButterAudio] = Butterworth(10, audioSamplingPeriod, audioN, audioTimeAxis, false);
audioFilteredSignal = filtfilt(bButterAudio, aButterAudio, audio);

figure;
plot(audioTimeAxis, audio);
xlabel('Tiempo [s]');
ylabel('Intensidad de audio');
print -dpng 'images/ej9/audio.png';

figure;

audioWindowWidth = ceil(5*1.2/audioSamplingPeriod);
audioSpecgramWindow = window(@hamming, audioWindowWidth);
generateSpectrograms(audio, audioFilteredSignal, audioSamplingPeriod, audioSpecgramWindow, audioWindowWidth);
print -dpng 'images/ej9/audioespectrograma.png';

audioDFT = fftshift(fft(audio));

% Calculo el eje de frecuencia discreta.
audioDiscreteFrequencyStep = 2*pi / audioN;
audioDiscreteFrequencyAxis = [-pi : audioDiscreteFrequencyStep : pi-audioDiscreteFrequencyStep];

% Calculo el eje de frecuencia angular y en Hertz.
audioContinuousFrequencyAxis = audioDiscreteFrequencyAxis / audioSamplingPeriod;
audioHertzAxis = audioContinuousFrequencyAxis / (2*pi);
plotDFT(audioHertzAxis, audioDFT);
print -dpng 'images/ej9/audiodft.png';





%% 10. Realizar un detector automático de latidos. El mismo debe tomar como entrada
% la señal FPG y producir como salida un vector de tiempos, donde cada tiempo
% corresponde a la detección de un latido en la señal. Para esto, se sugiere implementar
% los siguientes pasos:
printf('\nEjercicio 10 y 12:\n');


function newPeaks = applyRuleA (currentPeaks, newPeak, newPeakValue, signal, samplingPeriod)
	if length(currentPeaks) == 0
		newPeaks = [newPeak];
	else
		delta = newPeak - currentPeaks(length(currentPeaks));
		if delta*samplingPeriod > 0.400
			newPeaks = [currentPeaks, newPeak];
		else
			if newPeakValue > signal(currentPeaks(end))
				newPeaks = [currentPeaks, newPeak];
			else
				newPeaks = currentPeaks;
			end
		end
	end
endfunction


function newPeaks = applyRuleB(currentPeaks, newPeak, signal, samplingPeriod, currentUmbral)
	newPeaks = [currentPeaks];
	if length(currentPeaks) >= 2
		lastIBI = currentPeaks(length(currentPeaks)) - currentPeaks(length(currentPeaks)-1);
		newIBI = newPeak - currentPeaks(length(currentPeaks));
		if newIBI > 1.5*lastIBI
			morePeaks = findPulses(currentPeaks(length(currentPeaks)-1), newPeak, currentUmbral/2, signal, samplingPeriod);
			for peak = morePeaks
				if (peak - newPeaks(length(newPeaks))) * samplingPeriod > 0.360
					newPeaks = [newPeaks, peak];
				endif
			endfor
		end
	end
endfunction


% Reglas:
% 	extra) Desde que la señal sube más que el umbral hasta que vuelve a bajar, sólo 
%	puede haber un pico, y es el mayor de todos los picos de ese intervalo.
%
%	a. Si dos latidos se detectaron con una separación
%	temporal menor a 200ms, sobrevive sólo aquel que corresponda al pico de la
%	señal mayor entre ambos.
%
%	b. Si el IBI instantáneo aumenta repentinamente en al
%	menos 1.5 veces entre muestras consecutivas, puede haberse perdido la
%	detección de un latido. Cuando este sea el caso, realice una nueva detección
%	de picos dentro del intervalo correspondiente, utilizando un umbral la mitad
%	del nominal. Si de esta manera se halla un nuevo pico, distanciado al menos
%	360ms de la detección precedente, entonces clasificarlo como latido.
%
function pulseXs = findPulses (firstI, lastI, umbral, normalizedFilteredSignal, samplingPeriod)
	deltas = diff(normalizedFilteredSignal);
	minusUmbral = normalizedFilteredSignal .- umbral;
	signed = minusUmbral ./ abs(minusUmbral);
	signDeltas = diff(signed);

	returnValue = [];
	higherThanUmbral = false;
	lastPeakFound = firstI;
	lastPeakFoundValue = umbral;

	for i = [firstI+1 : lastI-1]
		
		% Máquina de estados
		if signDeltas(i) == 2
			% Superamos el umbral
			higherThanUmbral = true;

		elseif signDeltas(i) == -2
			% Cuando terminé de mirar un intervalo mayor que el umbral
			returnValue = applyRuleB(returnValue, lastPeakFound, normalizedFilteredSignal, samplingPeriod, umbral);
			returnValue = applyRuleA(returnValue, lastPeakFound, lastPeakFoundValue, normalizedFilteredSignal, samplingPeriod);
			lastPeakFoundValue = umbral;
			higherThanUmbral = false;
		endif
		
		% Estamos más alto que el umbral
		if higherThanUmbral

			% Estamos en un maximo local
			if (deltas(i-1) > 0) && (deltas(i+1) < 0)

				% Y es el más grande desde que pasamos el umbral
				if normalizedFilteredSignal(i) > lastPeakFoundValue
					lastPeakFoundValue = normalizedFilteredSignal(i);
					lastPeakFound = i;
				endif
			endif
		endif
	endfor
	pulseXs = returnValue;
endfunction


function outputSignal = myFilter (filterImpulseResponse, inputSignal)
	outputSignal = conv(filterImpulseResponse, inputSignal);
	outputSignal = outputSignal(floor(length(filterImpulseResponse)/2):end-ceil(length(filterImpulseResponse)/2));
endfunction


function signalPeaksTimes = detectPulses (signal, samplingPeriod, timeAxis)
	N = length(signal);

	% a. Filtrado pasa-banda de la señal, utilizando el filtrado del ejercicio 7.
	[bButter aButter] = Butterworth (10, samplingPeriod, N, timeAxis, false);
	iirFilteredSignal = filtfilt(bButter, aButter, signal);

	% b. Filtro de derivada, implementado con un filtro FIR h(n)=[-2 -1 0 1 2].
	hDriftRemoval = [-2 -1 0 1 2];
	filteredSignal = myFilter(hDriftRemoval, iirFilteredSignal);

	% c. Normalización con energía instantánea:

	% * Primero calcular la energía instantánea de la señal mediante un filtro MA 1
	% de la señal del punto 10a elevada al cuadrado
	maN = ceil(4/samplingPeriod);
	maFilter = ones(maN, 1) / maN;
	energy = myFilter(maFilter, iirFilteredSignal.^2);

	% luego dividir la señal del punto b por el vector obtenido.
	% Esto tiene como objeto reducir el impacto de la presión sanguínea sobre el
	% nivel de señal.
	% NOTA: Divido por la raiz cuadrada.
	filteredSignal = filteredSignal ./ sqrt(energy);

	% d. Sobre-muestreo en un factor 4 para obtener mayor resolución temporal:
	% implemente el sobre-muestreo utilizando la función upsample y diseñe un
	% filtro interpolador FIR utilizando la herramienta fdatool de Matlab.
	upsampledSignal = upsample(filteredSignal, 4);
	
	% Uso una sinc de 25 puntos
	hInterpolator = fir1(24, 1/4, 'low')/(1/4);

	interpolatedSignal = myFilter(hInterpolator, upsampledSignal);
	normalizedFilteredSignal = interpolatedSignal / max(interpolatedSignal);
	upsampledTimeAxis = [0:samplingPeriod/4:(samplingPeriod) * (N-1/4)];

	figure;
	plot(upsampledTimeAxis, interpolatedSignal);
	hold on;
	stem(timeAxis, filteredSignal);
	title('Interpolacion');
	legend('Señal interpolada', 'Señal original');
	xlabel('Tiempo [s]');
	ylabel('Intensidad');
	xlim([55 60]);
	ylim([-10 10]);
	print -dpng 'images/ej10/interpolacion.png';

	% Grafique: respuesta en frecuencia del filtro en módulo y fase, y señal original
	% y sobre-muestreada en superposición.
	figure;
	stem(timeAxis(1:length(hInterpolator)), hInterpolator);
	grid on;
	set(gca, 'xminorgrid', 'on');
	title('Respuesta al impulso');
	xlabel('Tiempo[s]');
	print -dpng 'images/ej10/impz.png';

	figure;
	freqz(hInterpolator, [1]);
	print -dpng 'images/ej10/freqz.png';

	% superposicion se plotea en el item siguiente

	% e. Detector de picos mediante umbral (puede definir como umbral un valor
	% arbitrario)
	% 12. Opcional. Mejorar el detector de latidos aplicando las reglas de [5]:
	% 	a. Establecer como regla que, si dos latidos se detectaron con una separación
	% 	temporal menor a 200ms, sobrevive sólo aquel que corresponda al pico de la
	% 	señal mayor entre ambos.
	%	edit: uso 400ms.
	% 	b. Establecer como regla que si el IBI instantáneo aumenta repentinamente en al
	% 	menos 1.5 veces entre muestras consecutivas, puede haberse perdido la
	% 	detección de un latido. Cuando este sea el caso, realice una nueva detección
	% 	de picos dentro del intervalo correspondiente, utilizando un umbral la mitad
	% 	del nominal. Si de esta manera se halla un nuevo pico, distanciado al menos
	% 	360ms de la detección precedente, entonces clasificarlo como latido.

	signalPeakXs = findPulses(1, 4*N, 0.015, normalizedFilteredSignal, samplingPeriod/4);
	signalPeaksTimes = upsampledTimeAxis(signalPeakXs);
	signalPeaksValues = normalizedFilteredSignal(signalPeakXs);

	figure;
	plot(timeAxis, signal);
	hold on;
	plot(upsampledTimeAxis, normalizedFilteredSignal*median(signal(ceil(signalPeakXs/4))));
	xlabel('Tiempo [s]');
	ylabel('Intensidad');
	grid on;
	legend('Señal', 'Señal filtrada');
	print -dpng 'images/ej10/superposicion.png';

	% f. Gráfico en superposición de la señal con las marcas de los picos detectados.
	figure;
	plot(upsampledTimeAxis, normalizedFilteredSignal);
	hold on;
	stem(signalPeaksTimes, signalPeaksValues);

	xlabel('Tiempo [s]');
	ylabel('Intensidad');
	ylim([-2 2]);
	xlim([64 74]);
	grid on;
	legend('Señal filtrada', 'Picos');
	print -dpng 'images/ej10/picos.png';

endfunction

pulseTimes = detectPulses(signal, samplingPeriod, timeAxis);





%% 11. En base a los resultados del punto anterior, calcule y grafique el intervalo
% temporal instantáneo entre latidos (IBI: inter-beat interval) y los LPM instantáneos.
printf('\nEjercicio 11:\n');

instantInterBreathIntervals = [0, diff(pulseTimes)];
figure;
stem(pulseTimes, instantInterBreathIntervals);
xlabel('Tiempo [s]');
ylabel('Inter breath interval [s]');
grid on;
print -dpng 'images/ej11/IBI.png';
