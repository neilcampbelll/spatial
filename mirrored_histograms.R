df <- tacsat %>%
  dplyr::select(LE_L5MET, SI_SP)

df2 <- df

df2$SI_SP <- -(df2$SI_SP)

df <- rbind(df, df2)

ggplot(df, aes(SI_SP)) + 
  geom_histogram(aes(fill = "VMS Pings"), alpha = 0.5, position = "identity", bins = 20) +
  facet_wrap( ~ LE_L5MET, scales = "free_y") +
  scale_fill_manual(values = c("VMS Pings" = "salmon")) +
  labs(title = "Histogram of Activities", x = "Speeds", y = "Frequency") +
  theme_minimal()
