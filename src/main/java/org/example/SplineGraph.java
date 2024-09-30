package org.example;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;

import static java.lang.Math.*;

public class SplineGraph extends JFrame {

    public SplineGraph(String title) {
        super(title);

        JPanel panel = new JPanel(new GridLayout(3, 1));

        JFreeChart tabulatedChart = createTabulatedChart();
        JFreeChart originalChart = createOriginalChart();
        JFreeChart splineChart = createSplineChart();
        JFreeChart errorChart = createErrorChart();

        ChartPanel tabulatedPanel = new ChartPanel(tabulatedChart);
        ChartPanel originalPanel = new ChartPanel(originalChart);
        ChartPanel splinePanel = new ChartPanel(splineChart);
        ChartPanel errorPanel = new ChartPanel(errorChart);

        panel.add(tabulatedPanel);
        panel.add(originalPanel);
        panel.add(splinePanel);
        panel.add(errorPanel);

        tabulatedPanel.setPreferredSize(new Dimension(800, 200));
        originalPanel.setPreferredSize(new Dimension(800, 200));
        splinePanel.setPreferredSize(new Dimension(800, 200));
        errorPanel.setPreferredSize(new Dimension(800, 200));

        setContentPane(panel);
    }

    private JFreeChart createTabulatedChart() {
        XYSeries tabulatedPoints = new XYSeries("Tabulated Points");

        int N = 50;
        double x0 = 0.0;
        double xn = 6.0;
        double hh = (xn - x0) / N;
        double[] x = new double[N + 1];
        double[] y = new double[N + 1];

        for (int i = 0; i <= N; i++) {
            x[i] = x0 + i * hh;
            y[i] = f(x[i]);
            tabulatedPoints.add(x[i], y[i]);
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(tabulatedPoints);

        return ChartFactory.createScatterPlot(
                "Tabulated Points",
                "X", "Y",
                dataset,
                PlotOrientation.VERTICAL,
                true, true, false);
    }

    private JFreeChart createOriginalChart() {
        XYSeries originalFunction = new XYSeries("Original Function f(x)");

        int N = 50;
        double x0 = 0.0;
        double xn = 6.0;
        double hh = (xn - x0) / N;

        for (int i = 0; i <= N; i++) {
            double xx = x0 + i * hh;
            originalFunction.add(xx, f(xx));
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(originalFunction);

        return ChartFactory.createXYLineChart(
                "Original Function f(x)",
                "X", "Y",
                dataset,
                PlotOrientation.VERTICAL,
                true, true, false);
    }

    private JFreeChart createSplineChart() {
        XYSeries splineFunction = new XYSeries("Spline Approximation");

        int N = 50;
        double x0 = 0.0;
        double xn = 6.0;
        double hh = (xn - x0) / N;
        double[] x = new double[N + 1];
        double[] y = new double[N + 1];
        double[] h = new double[N + 1];
        double[] a = new double[N + 1];
        double[] b = new double[N + 1];
        double[] c = new double[N + 1];
        double[] d = new double[N + 1];

        for (int i = 0; i <= N; i++) {
            x[i] = x0 + i * hh;
            y[i] = f(x[i]);
            h[i] = hh;
        }

        progonka(y, h, N, c);

        for (int i = 1; i < N; i++) {
            a[i] = y[i - 1];
            b[i] = (y[i] - y[i - 1]) / h[i] - (h[i] / 3) * (c[i + 1] + 2 * c[i]);
            d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
        }
        a[N] = y[N - 1];
        b[N] = (y[N] - y[N - 1]) / h[N] - (2.0 / 3.0) * h[N] * c[N];
        d[N] = -c[N] / (3 * h[N]);

        for (int i = 0; i <= N - 1; i++) {
            for (double xx = x[i]; xx < x[i + 1]; xx += 0.01) {
                double splineValue = a[i + 1] + b[i + 1] * (xx - x[i]) + c[i + 1] * pow((xx - x[i]), 2)
                        + d[i + 1] * pow((xx - x[i]), 3);
                splineFunction.add(xx, splineValue);
            }
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(splineFunction);

        saveToFile(x, y, h, a, b, c, d, N);

        return ChartFactory.createXYLineChart(
                "Spline Approximation",
                "X", "Y",
                dataset,
                PlotOrientation.VERTICAL,
                true, true, false);
    }

    public static void progonka(double[] y, double[] h, int N, double[] c) {
        double[] alfa = new double[N + 1];
        double[] beta = new double[N + 1];
        double[] hamma = new double[N + 1];
        double[] delta = new double[N + 1];
        double[] A = new double[N + 1];
        double[] B = new double[N + 1];

        alfa[1] = hamma[1] = delta[1] = 0.0;
        beta[1] = 1.0;

        for (int i = 2; i <= N; i++) {
            alfa[i] = h[i - 1];
            beta[i] = 2 * (h[i - 1] + h[i]);
            hamma[i] = h[i];
            delta[i] = 3 * (((y[i] - y[i - 1]) / h[i]) - ((y[i - 1] - y[i - 2]) / h[i - 1]));
        }
        hamma[N] = 0.0;

        A[1] = -hamma[1] / beta[1];
        B[1] = delta[1] / beta[1];
        for (int i = 2; i <= N - 1; i++) {
            A[i] = -hamma[i] / (alfa[i] * A[i - 1] + beta[i]);
            B[i] = (delta[i] - alfa[i] * B[i - 1]) / (alfa[i] * A[i - 1] + beta[i]);
        }

        c[N] = (delta[N] - alfa[N] * B[N - 1]) / (alfa[N] * A[N - 1] + beta[N]);
        for (int i = N; i > 1; i--) {
            c[i - 1] = A[i - 1] * c[i] + B[i - 1];
        }
    }


    public static void saveToFile(double[] x, double[] y, double[] h, double[] a, double[] b, double[] c, double[] d, int N) {
        try (PrintWriter writer = new PrintWriter(new FileWriter("output.txt"))) {
            double[] xm = new double[20 * N + 1];
            double[] ym = new double[20 * N + 1];

            double hh = h[0] / 20.0;
            for (int i = 0; i <= 20 * N; i++) {
                xm[i] = x[0] + i * hh;
                ym[i] = f(xm[i]);
            }

            double s, eps;
            int j = 1;
            for (int i = 0; i <= 20 * (N - 1); i++) {
                s = a[j] + b[j] * (xm[i] - x[j - 1]) + c[j] * pow((xm[i] - x[j - 1]), 2)
                        + d[j] * pow((xm[i] - x[j - 1]), 3);
                eps = abs(s - ym[i]);
                writer.printf("[%d,%d]\t%e\t%e\t%e\t%e%n", i, j, xm[i], ym[i], s, eps);
                if ((i != 0) && (i % 20 == 0)) {
                    j++;
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private JFreeChart createErrorChart() {
        XYSeries errorSeries = new XYSeries("Error of Spline Approximation");

        int N = 50;
        double x0 = 0.0;
        double xn = 6.0;
        double hh = (xn - x0) / N;
        double[] x = new double[N + 1];
        double[] y = new double[N + 1];
        double[] h = new double[N + 1];
        double[] a = new double[N + 1];
        double[] b = new double[N + 1];
        double[] c = new double[N + 1];
        double[] d = new double[N + 1];

        for (int i = 0; i <= N; i++) {
            x[i] = x0 + i * hh;
            y[i] = f(x[i]);
            h[i] = hh;
        }

        progonka(y, h, N, c);

        for (int i = 1; i < N; i++) {
            a[i] = y[i - 1];
            b[i] = (y[i] - y[i - 1]) / h[i] - (h[i] / 3) * (c[i + 1] + 2 * c[i]);
            d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
        }
        a[N] = y[N - 1];
        b[N] = (y[N] - y[N - 1]) / h[N] - (2.0 / 3.0) * h[N] * c[N];
        d[N] = -c[N] / (3 * h[N]);

        for (int i = 0; i <= N - 1; i++) {
            for (double xx = x[i]; xx < x[i + 1]; xx += 0.01) {
                double splineValue = a[i + 1] + b[i + 1] * (xx - x[i]) + c[i + 1] * pow((xx - x[i]), 2)
                        + d[i + 1] * pow((xx - x[i]), 3);
                double errorValue = Math.abs(splineValue - f(xx));
                errorSeries.add(xx, errorValue);
            }
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(errorSeries);

        return ChartFactory.createXYLineChart(
                "Error of Spline Approximation",
                "X", "Error",
                dataset,
                PlotOrientation.VERTICAL,
                true, true, false);
    }

    public static double f(double x) {
        return sin(x);
    }

    public static void main(String[] args) {
        SplineGraph chart = new SplineGraph("Spline Interpolation Example");
        chart.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        chart.pack();
        chart.setLocationRelativeTo(null);
        chart.setVisible(true);
    }
}
