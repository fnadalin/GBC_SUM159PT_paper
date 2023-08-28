
s_names <- c("top 1%", "top 10%", "top 50%", "top 100%")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
        cat("\nUsage: Rscript plot_cumulative_cpm.R <cpm.tsv> <case> <control> <out>\n")
        q()
}

COUNTS <- args[1]
CASE <- args[2]
CTRL <- args[3]
OUT <- args[4]

data <- read.table(COUNTS, sep = "\t", header = TRUE)

case <- data[[CASE]]
ctrl <- data[[CTRL]]
names(case) <- names(ctrl) <- rownames(data)
case <- case[ctrl > 0]

count_case <- sort(case, decreasing = TRUE)
cum <- count_case
for (i in 2:length(cum)) {
    cum[i] <- cum[i-1]+count_case[i]
}
n <- length(cum)
s <- c(floor(n/100),floor(n/10),floor(n/2),n)
x <- c(cum[s[1]], cum[s[2]]-cum[s[1]], cum[s[3]]-cum[s[2]], cum[s[4]]-cum[s[3]])
cum_case <- cum

count_ctrl <- ctrl[match(names(count_case),names(ctrl))]
cum <- count_ctrl
for (i in 2:length(cum)) {
    cum[i] <- cum[i-1]+count_ctrl[i]
}
y <- c(cum[s[1]], cum[s[2]]-cum[s[1]], cum[s[3]]-cum[s[2]], cum[s[4]]-cum[s[3]])
cum_ctrl <- cum

library("scales")
colors <- c("blue", alpha("blue", alpha = 0.6), alpha("blue", alpha = 0.4), alpha("blue", alpha = 0.2))

TITLE <- paste(CTRL, CASE, sep="-")
pdf(OUT, width = 6.5, height = 6.5)
par(fig = c(0.2,1,0.2,1))
plot(cum_ctrl,cum_case, axes = FALSE, xlab = NA, ylab = NA, type = "l", lwd = 4, main = TITLE)
legend("bottomright", legend = s_names, bty = "n", fill = colors)
par(fig = c(0.12,0.37,0.21,0.98), new = TRUE)
barplot(matrix(x,nrow=4,ncol=1), ylab = "cumulative CPM (case)", col = colors)
par(fig = c(0.215,0.98,0.07,0.41), new = TRUE)
barplot(matrix(y,nrow=4,ncol=1), horiz = TRUE, xlab = "cumulative CPM (control)", col = colors)
axis(3, at = cum_ctrl[s], labels = s)
dev.off()

q()


