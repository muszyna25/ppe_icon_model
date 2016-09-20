\LoadClass[11pt,a4paper,titlepage]{scrreprt}
\RequirePackage[english]{babel}
\RequirePackage{fancyhdr}
\RequirePackage{natbib}
%\RequirePackage{chicago}
\RequirePackage{fancybox}
\RequirePackage{fancyvrb}
\RequirePackage{colordvi}
\RequirePackage{supertabular}
\RequirePackage{longtable}
\RequirePackage{lscape}       %zur Verwendung von Querformat-Blöcken
\RequirePackage{multirow}
\RequirePackage{threeparttable}
\RequirePackage{randtext}
\RequirePackage{array}
\RequirePackage{ifpdf}
\RequirePackage{amsmath}
\RequirePackage{subfig}
\RequirePackage{float}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\ifpdf
\RequirePackage{hyperref}
\RequirePackage[pdftex]{graphicx,color,rotating}
\definecolor{UrlColor}{rgb}{0,0,0.95}
\definecolor{LinkColor}{rgb}{0.75,0,0}
\definecolor{CiteColor}{rgb}{0,0.7,0}
\definecolor{MenuColor}{rgb}{0.2,0.8,0}
\definecolor{FileColor}{rgb}{0,0.7,0.7}
% Set PDF document information (optional)
\hypersetup{%
  pdftitle={ICON User's Guide},
  pdfauthor={Max Mustermann},
  pdfsubject={PDFLaTeX ICON User's Guide},
  pdfkeywords={PDFLaTeX},
  colorlinks=true,linkcolor=LinkColor,urlcolor=UrlColor %
                 ,citecolor=CiteColor,menucolor=MenuColor %
                 ,filecolor=FileColor
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage[dvips]{graphicx,color,rotating}
\DeclareGraphicsExtensions{.eps,.eps2,.ps}
\def\hypertarget#1#2{\,{#2}\,}
\def\hyperlink#1#2{\,{#2}\,}
\def\href#1#2{\,{#2}\,}
\fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Configuration
% \setlength{\unitlength}{1mm}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\voffset-1.0cm
\hoffset-1.0cm
%
\setlength{\textwidth}{16.4cm}
\setlength{\topmargin}{0.8cm}
\setlength{\headheight}{0.5cm}
\setlength{\headsep}{0.6cm}
\setlength{\parindent}{0.0cm}
\setlength{\parskip}{0.3cm}
\setlength{\textheight}{23.7cm}
\setlength{\footskip}{1.4cm}
%
\renewcommand{\floatpagefraction}{0.95}
\renewcommand{\textfraction}{0.05}
\renewcommand{\topfraction}{0.9}
\renewcommand{\bottomfraction}{0.5}
%
\frenchspacing
%
\newenvironment{Ventry}[1]%
  {\begin{list}{}{\renewcommand{\makelabel}[1]{\quad\textsf{##1}\hfil}%
     \settowidth{\labelwidth}{\textsf{\quad#1\quad}}%
     \setlength{\leftmargin}{\labelwidth}}}%
  {\end{list}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enable use of external hyperlinks
\usepackage{url}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\graphicspath{{./}{pictures/}}

\newcolumntype{C}[1]{>{\centering\arraybackslash}p{#1}}

\usepackage[blocks]{authblk}
\usepackage{makeidx}
\makeindex
\pagenumbering{gobble}
\newcommand{\contrib}[1]{\textit{Authors of the original document:\\ {#1}}}
\newcommand{\prt}{\partial}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for code listings:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage{listings}
\definecolor{dgray}{gray}{0.9}
\lstset{basicstyle=\normalsize\ttfamily,
%         frame=lines,
        backgroundcolor=\color{dgray},
        breaklines=true,
        language=fortran,
        escapechar=?
        } 
