using Weave

reportLocation = "/Users/adamaker/Desktop/Research/Thesis/ThesisCode/Thesis/Report.jml"
publishLocation = "/Users/adamaker/Desktop/Research/Thesis/ThesisReports/Report.pdf"
weave(reportLocation, out_path = publishLocation, doctype = "pandoc2pdf")