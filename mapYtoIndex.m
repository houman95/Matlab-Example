function index = mapYtoIndex(Y)
    switch Y
        case '0'
            index = 1;
        case '1'
            index = 2;
        case 'C'
            index = 3;
        case 'I'
            index = 4;
        case '+'
            index = 5;
        case '-'
            index = 6;
        otherwise
            error('Invalid Y value');
    end
end
