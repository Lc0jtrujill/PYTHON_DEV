function y = v10_exact_smooth(t)
global V10_ALPHA V10_T
y = 1 + (t./V10_T).^(5+V10_ALPHA);
end
