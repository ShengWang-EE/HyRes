function dailyProfile = getDailyProfile(iDay,historicalYearData)
% 根据天数从全年数据中提取。这样只用从内存中提取，节约时间，可以放在for循环中
startHour = (iDay-1) * 24 + 1;
endHour = iDay * 24;

dailyProfile.electricityDemand = historicalYearData.electricityDemand(startHour:endHour);
dailyProfile.gasDemand = historicalYearData.gasDemand(startHour:endHour);
dailyProfile.renewableCapacity = historicalYearData.renewableCapacity(startHour:endHour);

end