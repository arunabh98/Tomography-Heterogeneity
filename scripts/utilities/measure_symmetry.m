function s = measure_symmetry(test_image)

    img_f = fftshift(fft(test_image));
    
    img_fr = real(img_f);
    img_fz = imag(img_f);
    
    s = 1 - norm(img_fr, 2).^2/(norm(img_fr, 2).^2 + norm(img_fz, 2).^2);
end