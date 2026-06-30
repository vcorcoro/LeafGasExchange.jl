# System for tracking time relative to circadian cycle
@system ZeitgeiberTime begin
    ztime(context.clock.time): zeitgeiber_time ~ track(u"hr")
    zt(ztime): zeitgeiber_time_hour => ztime % 24u"hr" ~ track(u"hr")
    zt0: start_of_circadian_cycle => 18 ~ preserve(parameter, u"hr")
    ctime(ztime, zt0): calendar_time => ztime + zt0 ~ track(u"hr")
    ct(ctime): calendar_time_hour => ctime % 24u"hr" ~ track(u"hr")
end