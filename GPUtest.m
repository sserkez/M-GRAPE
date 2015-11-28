%lfla;g

A1 = rand(5000,5000);
tic;
B1 = fft(A1);
time1 = toc;

A2 = gpuArray(A1);
tic;
B2 = fft(A2);
time2 = toc;

speedUp = time1/time2;
disp(speedUp)

%%

tic;
A4 = rand(3000,3000);
B4 = fft(A4);
time4 = toc;

tic;
A5 = parallel.gpu.GPUArray.rand(3000,3000);
B5 = fft(A5);
B5 = gather(B5);
time5 = toc;

speedUp = time4/time5;
disp(speedUp);