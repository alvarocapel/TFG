function numBusinessDays = countBusinessDaysColumn(startDates, endDates)
    % This function calculates the number of business days between two columns of dates
    % Excluding weekends (Saturday and Sunday)
    
    % Create an array of business days for each pair of start and end dates
    numBusinessDays = arrayfun(@(startDate, endDate) sum(weekday(startDate:endDate) ~= 1 & weekday(startDate:endDate) ~= 7), startDates, endDates);
end
