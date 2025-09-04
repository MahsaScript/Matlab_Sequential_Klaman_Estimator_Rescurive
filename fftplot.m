function [f,y]=fftplot(x,step)

[r,c] = size(x);

if r == 1
	x = x.';   % make it a column
end;

[n,cc] = size(x);

m = 2^nextpow2(n);
y = fft(x,m);
y = y(1:n,:);
yy=abs(y)*step;
kk=1:m/2+1;
f(kk)=(kk-1)/m/step;

 
y=yy(1:m/2+1);
plot(f,y);
xlabel('Frequency (Hz)');
ylabel('Fourier Amplitude');
%title('Fouier Transform of the Data')

