data_a <- read.table("exp_4_amplicon_v2/4A_DNA/D15_A/out_d1_f1/groups.tsv", sep = "\t", header = TRUE)
data_b <- read.table("exp_4_amplicon_v2/4A_DNA/D15_B/out_d1_f1/groups.tsv", sep = "\t", header = TRUE)
data_c <- read.table("exp_4_amplicon_v2/4A_DNA/D15_C/out_d1_f1/groups.tsv", sep = "\t", header = TRUE)

all_gbc <- unique(c(as.character(data_a$hub), as.character(data_b$hub), as.character(data_c$hub)))
df <- data.frame(gbc_id = all_gbc, a = rep(0, length(all_gbc)), b = rep(0, length(all_gbc)), c = rep(0, length(all_gbc)))
df$a[match(data_a$hub, all_gbc)] <- data_a$hub_count
df$b[match(data_b$hub, all_gbc)] <- data_b$hub_count
df$c[match(data_c$hub, all_gbc)] <- data_c$hub_count

write.table(df, file = "EXP4A_bulk_GBC.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

pdf("EXP4A_bulk_GBC.pdf", width = 9, height = 3)
par(mfrow=c(1,3))
plot(df$a,df$b,main="")
plot(df$a,df$c,main="EXP4A")
plot(df$b,df$c,main="")
dev.off()

data_a <- read.table("exp_1_amplicon_v2/1C_DNA/D15_A/out_d1_f1/groups.tsv", sep = "\t", header = TRUE)
data_b <- read.table("exp_1_amplicon_v2/1C_DNA/D15_B/out_d1_f1/groups.tsv", sep = "\t", header = TRUE)
data_c <- read.table("exp_1_amplicon_v2/1C_DNA/D15_C/out_d1_f1/groups.tsv", sep = "\t", header = TRUE)

all_gbc <- unique(c(as.character(data_a$hub), as.character(data_b$hub), as.character(data_c$hub)))
df <- data.frame(gbc_id = all_gbc, a = rep(0, length(all_gbc)), b = rep(0, length(all_gbc)), c = rep(0, length(all_gbc)))
df$a[match(data_a$hub, all_gbc)] <- data_a$hub_count
df$b[match(data_b$hub, all_gbc)] <- data_b$hub_count
df$c[match(data_c$hub, all_gbc)] <- data_c$hub_count

write.table(df, file = "EXP1C_bulk_GBC.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

pdf("EXP1C_bulk_GBC.pdf", width = 9, height = 3)
par(mfrow=c(1,3))
plot(df$a,df$b,main="")
plot(df$a,df$c,main="EXP1C")
plot(df$b,df$c,main="")
dev.off()

