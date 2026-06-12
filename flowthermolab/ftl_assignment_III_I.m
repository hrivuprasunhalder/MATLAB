numbers = zeros(1,10);

for i = 1:10
    numbers(i)  = input('Enter number :');
end    

sorted = sort(numbers);

disp(' Numbers in increasing order : ')
disp(sorted)