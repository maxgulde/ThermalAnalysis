% Convert STK date format to MATLAB date string
function date = getDate(stkDate)

% stkDate = '22 Jul 2018 13:48:58.676';

date = datetime(stkDate,'InputFormat','dd MMM yyyy HH:mm:ss.SSS');