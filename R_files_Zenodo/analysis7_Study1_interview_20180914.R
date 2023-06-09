library(ggplot2) # きれいなグラフをプロットするためのライブラリを読み込み
library(signal) # ハイカットフィルタのためのライブラリを読み込み
library(gridExtra) # グラフを縦横に並べるためのライブラリを読み込み
library(ggfortify) # 相互相関関数をプロットして並べるためのライブラリを読み込み
library(gplots) # ヒートマップを表示するためのライブラリを読み込み
library(reshape2) #
library(pagenum) # ページ番号をつける

sampling.rate <- 7.5 # サンプリング周波数は7.5 Hz
nyquist.freq <- sampling.rate/2 # ナイキスト周波数は定義からサンプリング周波数の半分
filter.order <- 4 # ハイカットに使うButterworthフィルタは4次のものを使う
high.cut.threshold <- 0.1 # 0.1 Hz以上を遮断
low.cut.threshold <- 0.01 # 0.01 Hz以下を遮断
ccf.lag.max <- 20*sampling.rate # 長い相互相関の表示は20秒まで
thetheme <- theme_classic(6) # この数字を変えるとフォントサイズが変わる
thepointsize <- 0.1 # この数字を変えると散布図の点のサイズが変わる
maxfignum <- 10 # 一ページあたりの図の最大数
onlyheatmap <- TRUE # ヒートマップ以外の散布図やcross-correlogramを表示しないならTRUE、表示するならFALSE
no.p.value <- FALSE # p値のCSVファイルを出力しないならTRUE
analyze.raw <- FALSE # rawデータもヒートマップを作って検定するならTRUE
np <- 31 # ヒストグラムの階調数
palette <- colorRampPalette(c("blue", "white", "red"))(n=np)
unified.y.scale <- FALSE # 時系列の縦軸の範囲をすべて同一にするならばTRUE、チャネルごとに適宜決めるならFALSE
plus.low.cut <- TRUE # ローカットも入れるならTRUE、入れないならFALSE
separated.r.t <- TRUE # 最後のグラフでrestとtaskを分けるならTRUE、分けないならFALSE
#offsetcount <- round(100/0.13) # taskの始めと終わりから除去する長さ
offsetcount <- 0 # taskの始めと終わりから除去する長さ
article.mode <- TRUE # 論文の図を生成する

setPagenum(1)
addpagenum <- function() {
    if (!article.mode) {
        pagenum(col="gray25", just=c("center", "bottom"))
    }
}

channels <- 40
load.and.split <- function(filename) {
    rawdata <- as.data.frame(read.table(filename, sep=",", dec=".", stringsAsFactors=F))
    offset <- 6
    elinenum <- which(rawdata[[1]]=='Event')
    originalch <- rawdata[elinenum-1,offset:(channels+offset-1)]
    rawdata <- rawdata[(elinenum+1):nrow(rawdata),]
    rawdata <- cbind(rawdata[1], as.data.frame(lapply(rawdata[2], as.numeric)), as.data.frame(lapply(rawdata[offset:(offset+channels-1)], as.numeric)))
    channelname <- gsub("[^a-z0-9]*", "", originalch)
    colnames(rawdata) <- c(c('Event', 'Time'), channelname)
    leftside <- 1:10
    rightside <- c(20,18,19,17,15,16,14,12,13,11)
    rawdata <- rawdata[c(1,2,leftside+2,rightside+2,leftside+20+2,rightside+20+2)]
    return(rawdata)
}
#files <- list('ID016.csv', 'ID017.csv')
#files <- list('ID016.csv')
files <- list('ID016.csv', 'ID017.csv', 'ID018.csv', 'ID019.csv', 'ID020.csv')
#files <- list('ID016.csv', 'ID017.csv', 'ID018.csv')

pdf("results.pdf") # PDFファイルに出力

rawcontents <- lapply(files, load.and.split) # ファイル名のリストfilesの中のファイルをすべて読み込んでcontentsに代入

names(rawcontents) <- files # contentsに読み込んだデータフレームそれぞれにファイル名のラベルをつける

cat("Filtering...\\n")

prefix <- function(x) { # ファイル名から末尾の.csvを消去し、-は_に変更
    return(gsub("-", "_", substr(x, 1, nchar(x)-4)))
}

getright <- function(x) {
    return(lapply(x, function(x) { return(strsplit(x, "\\\\.")[[1]][2]) }))
}

getleft <- function(x) {
    return(lapply(x, function(x) { return(strsplit(x, "\\\\.")[[1]][1]) }))
}

getnum <- function(x) {
    return(gsub("[^0-9]", "", x))
}

bfilter1 <- butter(filter.order, high.cut.threshold/nyquist.freq, type="low") # Butterworthフィルタを作成
if (plus.low.cut) {
    bfilter2 <- butter(filter.order, low.cut.threshold/nyquist.freq, type="high") # Butterworthフィルタを作成
    ff <- function(x){return(filtfilt(bfilter2,filtfilt(bfilter1,x)))}
} else {
    ff <- function(x){return(filtfilt(bfilter1,x))}
}
filtered.rawcontents <- lapply(rawcontents, function(data){return(data.frame(cbind(data[1:2], lapply(data[3:(channels+2)], ff))))})

lrstrip <- function(x) {
    return(sub("[[:blank:]]*$", "", sub("^[[:blank:]]*", "", x)))
}

contents <- list()
filtered.contents <- list()
thetime <- list()
name.of.contents <- list()
for (file in files) {
    theraw <- rawcontents[[file]]
    thefiltered <- filtered.rawcontents[[file]]
    nstart <- which(sapply(theraw[1], lrstrip)=="Noise-start")
    nend <- which(sapply(theraw[1], lrstrip)=="Noise-end")
    if (length(nstart)!=length(nend)) {
        cat("Noise-start and Noise-end do not balance.\\n")
        stopifnot(FALSE)
    }
    if (length(nstart)>0) {
        nsn <- cbind(nstart,nend)
        skip <- 0
        for (i in 1:nrow(nsn)) {
            theraw <- theraw[c(1:(nsn[i,1]-skip-1), (nsn[i,2]-skip+1):nrow(theraw)),]
            skip <- skip-(nsn[i,2]-nsn[i,1]+1)
        }
    }
    rstart <- which(sapply(theraw[1], lrstrip)=="Rest-start")
    rend <- which(sapply(theraw[1], lrstrip)=="Rest-end")
    tstart <- which(sapply(theraw[1], lrstrip)=="Task-start")+offsetcount
    tend <- which(sapply(theraw[1], lrstrip)=="Task-end")-offsetcount
    istart <- which(sapply(theraw[1], lrstrip)=="Interview-start")+offsetcount
    iend <- which(sapply(theraw[1], lrstrip)=="Interview-end")-offsetcount
    cut.and.push <- function(start, end, r1, r2, clth, lr) {
        ind <- length(contents)+1
        contents[[ind]] <<- data.frame(theraw[start:end,(r1+2):(r2+2)])
        filtered.contents[[ind]] <<- data.frame(thefiltered[start:end,(r1+2):(r2+2)])
        thetime[[ind]] <<- theraw[start:end,2]
        name.of.contents[ind] <<- paste(clth, getnum(file), lr, sep="")
    }
    cut.and.push(rstart, rend, 1, 10, "rCl", "L")
    cut.and.push(rstart, rend, 11, 20, "rCl", "R")
    cut.and.push(tstart, tend, 1, 10, "tCl", "L")
    cut.and.push(tstart, tend, 11, 20, "tCl", "R")
    cut.and.push(istart, iend, 1, 10, "iCl", "L")
    cut.and.push(istart, iend, 11, 20, "iCl", "R")

    cut.and.push(rstart, rend, 21, 30, "rTh", "L")
    cut.and.push(rstart, rend, 31, 40, "rTh", "R")
    cut.and.push(tstart, tend, 21, 30, "tTh", "L")
    cut.and.push(tstart, tend, 31, 40, "tTh", "R")
    cut.and.push(istart, iend, 21, 30, "iTh", "L")
    cut.and.push(istart, iend, 31, 40, "iTh", "R")
}

names(contents) <- name.of.contents
names(filtered.contents) <- name.of.contents
names(thetime) <- name.of.contents

if (unified.y.scale) {
    ymax <- max(as.data.frame(lapply(contents, max))) # 読み込んだ値の最大値
    ymin <- min(as.data.frame(lapply(contents, min))) # 読み込んだ値の最小値
    yrange = ylim(ymin, ymax)
} else {
    yrange = c()
}

for (name1 in name.of.contents) {
    df <- cbind(thetime[[name1]], filtered.contents[[name1]])
    colnames(df) <- c("Time", colnames(filtered.contents[[name1]]))
    write.csv(df, paste("filtered_timeseries_", name1, ".csv", sep=""), row.names=FALSE)
}

the.name <- function(x) {
    return(substr(x, 2, nchar(x)-1))
}

is.self.simul <- function(x, y) {
    return(substr(x, 1, nchar(x)-1)==substr(y, 1, nchar(y)-1))
}

is.simul <- function(x, y) {
    return((substr(x, 1, 1)==substr(y, 1, 1)) && (gsub("[^0-9]*", "", x)==gsub("[^0-9]*", "", y)))
}

hemisphere <- function(x) {
    return(substr(x, nchar(x), nchar(x)))
}

trstate <- function(x) {
    return(substr(x, 1, 1))
}

fillzero <- function(x) {
    if (nchar(x)==1) {
        return(paste("0", x, sep=""))
    }
    return(x)
}

classify <- function(p) {
    p2 <- strsplit(p, "_")[[1]]
    x <- p2[[1]]
    y <- p2[[2]]
    if (is.self.simul(x, y)) {
        return("Intra-brain")
    } else if (is.simul(x, y)) {
        return("Inter-brain")
    }
    return("Control")
}

tile <- function(thelist, n1, n2, maxfignum) {
    N1 = ceiling(n1/maxfignum)
    N2 = ceiling(n2/maxfignum)
    for (i1 in 1:N1) {
        for (i2 in 1:N2) {
            newlist <- list()
            max1 <- 0
            max2 <- 0
            for (j1 in 1:maxfignum) {
                for (j2 in 1:maxfignum) {
                    pos1 <- (i1-1)*maxfignum+j1-1
                    pos2 <- (i2-1)*maxfignum+j2-1
                    if ((pos1<n1) && (pos2<n2)) {
                        newlist <- c(newlist, list(thelist[[pos1*n2+pos2+1]]))
                        max1 <- max(max1, j1)
                        max2 <- max(max2, j2 )
                    }
                }
            }
            do.call(grid.arrange, c(newlist, list(ncol=max2, nrow=max1))) # 散布図を縦max1行、横max2列の中に表示
        }
    }
}

add.new.row <- function(x, y)
{
    if (x==0) {
        return(y)
    }
    return(data.frame(x, y))
}

NG.channel.pair <- function(col1, col2) {
    l2r <- function(ch) {
        if (ch>10) {
            return(c(10,8,9,7,5,6,4,2,3,1)[ch-10])
        }
        return(ch)
    }
    chnum1 <- l2r((as.numeric(tail(strsplit(col1, "h")[[1]], n=1))-1)%%20+1)
    chnum2 <- l2r((as.numeric(tail(strsplit(col2, "h")[[1]], n=1))-1)%%20+1)
    checker <- rbind(c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
                     c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
                     c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
                     c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
                     c(FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE),
                     c(FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE),
                     c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
                     c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE),
                     c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE),
                     c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
    return(checker[chnum1,chnum2])
}

for (thecontents in list(list("raw", contents), list("filtered", filtered.contents))) { # 以下の内容をxcontents=contentsとxcontents=filtered.contentsでそれぞれ実行
    mode <- thecontents[[1]]
    xcontents <- thecontents[[2]]
    for (name1 in name.of.contents) { # 以下の内容をfilesの中のすべてのファイルについて実行
        cat(sprintf("Processing %s %s...\\n", mode, name1))
        dat1 <- xcontents[[name1]] # xcontentsの中に入っているファイル名がfile1と一致するデータフレームをdat1に代入
        n1 <- ncol(dat1) # このデータフレームのコラム数
        timeseries.plots <- list() # 時系列グラフを入れておくリストを初期化
        for (col1 in colnames(dat1)) {
            timeseries.plots <- c(timeseries.plots, list(ggplot(data=transform(dat1, sec=thetime[[name1]]), aes_string(x='sec',y=col1))+geom_line()+thetheme+yrange+labs(title=paste(name1, mode, sep=" ")))) # リストtimeseries.plotにグラフを追加する
	    if (length(timeseries.plots)==maxfignum) {
	   # do.call(grid.arrange, c(timeseries.plots, list(ncol=1, nrow=maxfignum))) # timeseries.plotsに入っているグラフを横1列、縦n1行に並べて表示
	    timeseries.plots <- list()
            addpagenum()
	    }
            # ch1とch6以外にsec（記録開始からの秒数）を追加したデータを使って描画する
            # 横軸がsec、縦軸がcol1
            # 折れ線グラフ
            # グラフはクラシックな体裁
            # 縦軸の最大値・最小値はymaxとymin
            # ファイル名をタイトルにする
        }
	if (length(timeseries.plots)>0) {
           # do.call(grid.arrange, c(timeseries.plots, list(ncol=1, nrow=maxfignum))) # timeseries.plotsに入っているグラフを横1列、縦n1行に並べて表示
            addpagenum()
}
    }
    if (analyze.raw | mode!="raw") {
    crosscorr.data.frame.self <- 0
    crosscorr.data.frame.simul <- 0
    crosscorr.data.frame.control <- 0
    r.corr <- c()
    t.corr <- c()
    i.corr <- c()
    for (name1 in name.of.contents) { # 以下の内容をfilesの中のすべてのファイルについて実行
        dat1 <- xcontents[[name1]] # xcontentsの中に入っているファイル名がfile1と一致するデータフレームをdat1に代入
        n1 <- ncol(dat1) # このデータフレームのコラム数（今回はch1とch6なのでn1=2になる）
        for (name2 in name.of.contents) { # 以下の内容をfilesの中のすべてのファイルについて実行
            if (match(name1, name.of.contents)<=match(name2, name.of.contents) && (trstate(name1)==trstate(name2) || !separated.r.t)) {
            cat(sprintf("Processing the pair of %s and %s...\\n", name1, name2))
            correlation.values = c()
            correlation.p = c()
            dat2 <- xcontents[[name2]] # xcontentsの中に入っているファイル名がfile1と一致するデータフレームをdat1に代入
            # 今回はfile1はid18ch1_6.csvのみなので、単にid18ch1_6.csvを読み込んだデータフレーム
            n2 <- ncol(dat2) # このデータフレームのコラム数（今回はch1とch6なのでn2=2になる）
            scatter.plots <- list() # 散布図のリストを初期化
            ccf.plots1 <- list() # 短い相互相関のリストを初期化
            ccf.plots2 <- list() # 長い相互相関のリストを初期化
            if (nrow(dat1)<nrow(dat2)) { # 行数の短い方に合わせて長い方を切る
                merged.data <- data.frame(dat1, head(dat2, n=nrow(dat1))) # 長さを合わせてマージ
            } else {
                merged.data <- data.frame(head(dat1, n=nrow(dat2)), dat2) # 長さを合わせてマージ
            }
            col1names <- lapply(colnames(dat1), function(x) { return(paste(name1, x, sep=".")) }) # ファイル名とチャンネル名をつなげた名前のリスト
            col2names <- lapply(colnames(dat2), function(x) { return(paste(name2, x, sep=".")) }) # ファイル名とチャンネル名をつなげた名前のリスト
            names(merged.data) <- c(col1names, col2names) # マージしたデータフレームに名前のリストをくっつける
            for (col1name in col1names) { # col1にch1とch6を代入して繰り返し
                for (col2name in col2names) { # col2にch1とch6を代入して繰り返し
                    #cat(sprintf("  %s and %s...\\n", col1name, col2name))
                    newplot <- ggplot(data=merged.data, aes_string(x=col1name,y=col2name))+geom_point(size=thepointsize)+thetheme
                    # 横軸がcol1、縦軸がcol2の散布図を作成
                    if ((name1==name2) && NG.channel.pair(col1name, col2name)) {
                        correlation.values <- c(correlation.values, NA)
                        correlation.p <- c(correlation.p, NA)
                    } else { # もしcol1とcol2が別物なら
                        result <- cor.test(merged.data[[col1name]], merged.data[[col2name]], method="spearman") # スピアマンの順位相関検定を実行
                        correlation.values <- c(correlation.values, result$estimate)
                        correlation.p <- c(correlation.p, result$p.value)
                        if (!onlyheatmap) {
                            newplot <- newplot+labs(title=sprintf("Spearman=%1.3f\\np=%1.3e", result$estimate, result$p.value)) # 相関係数とp値をタイトルとして表示（%1.3fは小数点より上が一桁、下が三桁の小数表示、%1.3eは小数点より上が一桁、下が3桁の指数表示）
                        }
                    }
                    if (!onlyheatmap) {
                        scatter.plots <- c(scatter.plots, list(newplot)) # 作った散布図を散布図のリストscatter.plotsに追加
                        ts1 <- ts(data=merged.data[[col1name]], frequency=sampling.rate) # col1のチャネルから時系列を作成
                        ts2 <- ts(data=merged.data[[col2name]], frequency=sampling.rate) # col2のチャネルから時系列を作成
                        thetitle <- sprintf("%s\\n%s", col1name, col2name) # グラフのタイトルを作成
                        ccf.plots1 <- c(ccf.plots1, list(autoplot(ccf(ts1, ts2, plot=FALSE))+thetheme+labs(title=thetitle))) # ラグの最大値をデフォルトにした場合の相互相関の図を作成するがこの場では表示せずccf.plots1に追加する
                        ccf.plots2 <- c(ccf.plots2, list(autoplot(ccf(ts1, ts2, lag.max=ccf.lag.max, plot=FALSE))+thetheme+labs(title=thetitle))) # ラグの最大値をccf.lag.maxとした場合の相互相関の図を作成するがこの場では表示せずccf.plots2に追加する
                    }
                }
            }
            if (!onlyheatmap) {
                cat(sprintf("  generating scatter plots...\\n"))
              #  tile(scatter.plots, n1, n2, maxfignum)
                cat(sprintf("  generating cross-correlogram A...\\n"))
              #  tile(ccf.plots1, n1, n2, maxfignum)
                cat(sprintf("  generating cross-correlogram B...\\n"))
               # tile(ccf.plots2, n1, n2, maxfignum)
            }
            correlation.values <- matrix(correlation.values, nrow=n2, ncol=n1)
            colnames(correlation.values) <- getright(col2names)
            rownames(correlation.values) <- getright(col1names)
            if (name1!=name2) {
                newrow <- as.data.frame(diag(correlation.values))
                names(newrow) <- paste(name1, "_", name2, sep="")
                if (is.self.simul(name1, name2)) {
                    crosscorr.data.frame.self <- add.new.row(crosscorr.data.frame.self, newrow)
                } else if (is.simul(name1, name2)) {
                    if (crosscorr.data.frame.simul==0) {
                        crosscorr.data.frame.simul <- newrow
                    } else {
                        crosscorr.data.frame.simul <- data.frame(crosscorr.data.frame.simul, newrow)
                    }
                } else {
                    if (crosscorr.data.frame.control==0) {
                        crosscorr.data.frame.control <- newrow
                    } else {
                        crosscorr.data.frame.control <- data.frame(crosscorr.data.frame.control, newrow)
                    }
                }
            }
            if (is.self.simul(name1, name2)) {
                c.val <- correlation.values
                if (name1==name2) {
                    c.val[lower.tri(c.val)] <- NA
                }
                if (hemisphere(name1)==hemisphere(name2) || (hemisphere(name1)=='L' && hemisphere(name2)=='R')) {
                    if (substr(name1, 1, 1)=="r") {
                        r.corr <- c(r.corr, na.omit(as.vector(c.val)))
                    } else if (substr(name1, 1, 1)=="t") {
                        t.corr <- c(t.corr, na.omit(as.vector(c.val)))
                    } else {
                        i.corr <- c(i.corr, na.omit(as.vector(c.val)))
                    }
                }
            }
            par(omi=c(3,0,0,3))
          #  heatmap.2(correlation.values, Rowv=FALSE, Colv=FALSE, key=TRUE, trace="none", denscol="black", na.color="black", col=palette, breaks=seq(-1,1,length=np+1), key.title=NA, density.info="none", dendrogram="none", key.xlab="Correlation", key.ylab="Count", key.par=par(cex=1, cex.axis=1, cex.label=1, mar=c(4,0,12,0)), xlab=substring(getleft(col2names), 2, 3), ylab=substring(getleft(col1names), 2, 3), main=paste(getleft(col2names)[[1]], getleft(col1names)[[1]]))
            addpagenum()
            
            correlation.p <- matrix(correlation.p, nrow=n2, ncol=n1)
            colnames(correlation.p) <- col2names
            rownames(correlation.p) <- col1names
            write.csv(correlation.values, paste("spearman_", mode, "_", name1, "_", name2, ".csv", sep=""))
            if (!no.p.value) {
                write.csv(correlation.p, paste("p_", mode, "_", name1, "_", name2, ".csv", sep=""))
            }
            }
        }
        if (length(r.corr)>0 && length(r.corr)==length(t.corr) && length(r.corr)==length(i.corr)) {
            x = 1:length(r.corr)
            lines <- cbind(r.corr, t.corr, i.corr)
            lcolor <- c("#ff2800", "#35a16b", "#0041ff")
            par(mar=c(3,3,3,6))
            matplot(x, lines, type="b", col=lcolor, pch=c(0,0,0), lwd=2, lty=1)
            par(xpd=TRUE)
            legend(par()$usr[2], par()$usr[4], legend=c("Rest", "Task", "Interview"), col=lcolor, pch=c(0,0,0), lwd=2, lty=1)
            test.rt <- wilcox.test(r.corr, t.corr, paired=T, alternative="t", exact=T)
            test.ri <- wilcox.test(r.corr, i.corr, paired=T, alternative="t", exact=T)
            test.ti <- wilcox.test(t.corr, i.corr, paired=T, alternative="t", exact=T)
            title(sprintf("%s pRT=%e pRI=%e pTI=%e", the.name(name1), test.rt$p.value, test.ri$p.value, test.ti$p.value), cex.main=0.5)
            addpagenum()
            r.corr <- c()
            t.corr <- c()
            i.corr <- c()
        }
    }
    chlabel <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10")

    if (separated.r.t) {
        lls <- c("r", "t", "i")
    } else {
        lls <- c("")
    }
    for (llabel in lls) {
        df <- data.frame(ch=chlabel, crosscorr.data.frame.self, crosscorr.data.frame.simul, crosscorr.data.frame.control)
        df <- melt(df)
        names(df) <- c("ch", "pairs", "correlation")
        df1 <- data.frame(stype=sapply(as.character(df$pairs), classify))
        df <- data.frame(df, df1)
        if (llabel!="") {
            df <- df[sapply(df$pairs, function(x) {return(length(strsplit(as.character(x), llabel)[[1]])==3)}),]
        }
        for (i in chlabel) {
            control <- subset(df, (ch==i) & (stype=="Control"), c(correlation))$correlation
            simultaneous <- subset(df, (ch==i) & (stype=="Inter-brain"), c(correlation))$correlation
            wresult <- wilcox.test(control, simultaneous, paired=FALSE, exact=TRUE)
            print(sprintf("%s ch%s, Cliff's delta=%e", llabel, i, 2*wresult$statistic/(length(control)*length(simultaneous))-1), quote=FALSE)
            df2 <- subset(df, ch==i, c(ch, pairs, correlation, stype))
            df2 <- transform(df2, stype=factor(stype, levels=c("Intra-brain", "Inter-brain", "Control")))
            if (article.mode) {
                g <- ggplot(df2, aes(x=stype, y=correlation, fill=pairs)) + geom_bar(position="dodge", stat="identity", show.legend=TRUE)+labs(title="")+labs(x="")+labs(y="Correlation coefficient")+labs(title=sprintf("%s ch%s, p=%e", llabel, i, wresult$p.value))+theme_grey(base_size=18)+theme(axis.text=element_text(size=24))+theme(legend.position="none")
            } else {
                g <- ggplot(df2, aes(x=stype, y=correlation, fill=pairs)) + geom_bar(position="dodge", stat="identity", show.legend=FALSE)+theme_grey(base_size=18)+theme(axis.text=element_text(size=18))+labs(title=sprintf("%s ch%s, p=%e", llabel, i, wresult$p.value))+xlab("")+ylab("Correlation coefficient")
            }
        plot(g)
        addpagenum()
    }
    }

    
    }
}

dev.off() # PDFとしてすべて出力
